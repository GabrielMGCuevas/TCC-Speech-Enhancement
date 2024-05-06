function [x,nurr,out,nur]=IMMIT(n_ipt,fs,win_t,ord,o2,iter)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo


%% Inicialização do filtro
len = fix(win_t * fs);

%% Cortando a entrada nas janelas que vamos trabalhar
len = fix(win_t * fs);
len_ovrlp=round(len/2);

y=buffer(n_ipt,len,len_ovrlp);
n=size(y);
delay=ord-1;

%% Definição de parêmtros do filtro
NumModel=3;
C{1} = [zeros(1,ord-1),1];
R{1} = o2;
P{1} = o2 * eye(ord);

x = zeros(1,length(n_ipt));
xa = zeros(ord,1);
xa(end)= n_ipt(1);

MarkovModel
nu=[0.5,0.25,0.25];
xo{1}=xa;

for ii=2:NumModel
    P{ii}=P{1};
    xo{ii}=xo{1};
    C{ii}=C{1};
    R{ii}=R{1};
end

xt=xo{1};
Pt=P{1};
x2=xo{1};
P2=P{1};

j = 1;

    %   Obtenção do modelo da primeira janela
    [a{2},o2w{2}] = KF_LPC(y(:,1),o2,ord,iter-1,x2,P2);
    a{3}=a{2};o2w{3}=o2w{2};
    

%% Processamento do sinal
for k = 1:n(2)

    %   Realocação dos modelos obtidos em janelas passadas
    a{1}=a{2};o2w{1}=o2w{2};
    a{2}=a{3};o2w{2}=o2w{3};
    %   Obtenção do modelo da janela futura por meio de um KF
    [a{3},o2w{3}] = KF_LPC2(y(:,k+1*(k<n(2))),o2,ord,iter-1,x2,P2);
    for ii=1:NumModel
        A{ii} = [zeros(ord-1,1), eye(ord-1); -flip(a{ii}(2:end))];
    end

    %...........................................................................................IMM
    for i = 1:len
        %% Reinicialização
        mu_=(Pr.*nu./(nu*Pr)')';

        for jj=1:NumModel
            sm=[zeros(size(xo{jj}))];
            for ii=1:NumModel
                sm=sm+mu_(ii,jj)*xo{ii};
            end
            xa_{jj}=sm;
        end
        for jj=1:NumModel
            sm=[zeros(size(P{ii}))];
            for ii=1:NumModel
                sm=sm + mu_(ii,jj)* (P{ii} +...
                    (xo{ii} - xa_{jj}) * (xo{ii} - xa_{jj})');
            end
            Pa_{jj}=sm;
        end
        clear sm

        %% Filtros de Kalman
        for m=1:NumModel
            x_{m} = A{m} * xa_{m};

            Pc{m} = (A{m} * Pa_{m} * A{m}') + ( o2w{m} * diag(C{m}));

            z{m} = y(i,k) - (C{m}*x_{m});

            S{m} = C{m}*Pc{m}*C{m}' + R{m};

            K{m} = (Pc{m} * C{m}')/S{m};

            xo{m} = x_{m} + K{m} * (z{m});

            P{m} = (eye(ord) - K{m} * C{m})* Pc{m};
        end

        %% Cálculo de Verossimilhança
        lkhd = zeros(NumModel,1);
        for m = 1:NumModel
            lkhd(m) = mvnpdf(z{m}, 0, S{m});
        end
        lkhd=lkhd+1e-50;
        nu = (Pr'*nu').*lkhd/( sum((Pr'*nu').*lkhd));
        nu = nu';
        %% Estimativa/Covariância Global
        %estimativa total de x
        xt=sum( (nu.*cell2mat(xo))' )';
        %estimativa total de P
        sm=[zeros(size(Pc{m}))];
        for m=1:NumModel
            sm=sm+nu(m) * (P{m} + (xo{m}-xt)*(xo{m}-xt)');
        end
        Pt=sm;
        clear sm

        xx(i,k)=xt(end);
        out(i,k)=xt(end-delay); %  Estimativa da saída com delay

        %   Armazenamento de variáveis para transição para janela futura
        if i==len-len_ovrlp
            xreinit{1} = xo{2};
            xreinit{2} = xo{3};
            xreinit{3} = xt;
            Preinit{1} = P{2};
            Preinit{2} = P{3};
            Preinit{3} = Pt;
            mureinit = nu;
        end

        %obtenção de gráfico de nu
        nur(i,k,1)=nu(1);
        nur(i,k,2)=nu(2);
        nur(i,k,3)=nu(3);
        j = j+1;

    end
    %...........................................................................................IMM
    x2=xt;
    P2=Pt;
    nu=[mureinit(2:3)/sum(mureinit(2:3)),0];
    xo=xreinit;
    P=Preinit;
end


%% Adição de janelas sobrepostas 
OVERLAP_ADD
x=speech(1+delay:numel(t)+delay);

%Obtenção de nu para gráfico
o=nur(:,:,1);
nurr(:,1)=o(1:end);
o=nur(:,:,2);
nurr(:,2)=o(1:end);
o=nur(:,:,3);
nurr(:,3)=o(1:end);

end