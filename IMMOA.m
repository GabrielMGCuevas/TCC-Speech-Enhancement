function [x,nurr]=IMMOA_Nghb(n_ipt,fs,win_t,ord_a,ord_b,o2,noisesample,iter,dataC)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo


%% Cortando a entrada nas janelas que vamos trabalhar
len = fix(win_t * fs);
len = fix(win_t * fs);
len_ovrlp=round(len/2);

y=buffer(n_ipt,len,len_ovrlp);
yc=buffer(dataC,len,len_ovrlp);
n=size(y);
delay=ord_a-1;


%% Inicialização do filtro
NumModel=3;
wndw=hamming(len);
C_a = [zeros(1,ord_a-1),1];
C_b= [zeros(1,ord_b-1),1];
C{1} = [zeros(1,ord_a-1),1,zeros(1,ord_b-1),1];
R{1} = o2;
P{1} = R{1} * eye(ord_a+ord_b);

x = zeros(1,length(n_ipt));
xa = zeros(ord_a+ord_b,1);

j = 1;

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

[a,o2w] = lpc(yc(:,1).*wndw, ord_a); %% Modelo de fala
[b,o2v] = lpc((y(:,1)-yc(:,1)).*wndw, ord_b); %% Modelo do ruído
A  = [zeros(ord_a-1,1), eye(ord_a-1);-flip(a(2:end))];
B = [zeros(ord_b-1,1), eye(ord_b-1);-flip(b(2:end))];
        
G{3}  = [A, zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), B];
        
Q{3}  = [o2w*diag(C_a), zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), o2v*diag(C_b)];

G{2}=G{3};Q{2}=Q{3};
%% Processamento do sinal
for k = 1:n(2)

    G{1}=G{2};Q{1}=Q{2};
    G{2}=G{3};Q{2}=Q{3};

    [a,o2w] = lpc(yc(:,k+1*(k<n(2))).*wndw, ord_a);
    [b,o2v] = lpc((y(:,k+1*(k<n(2)))-yc(:,k+1*(k<n(2)))).*wndw, ord_b);

    A  = [zeros(ord_a-1,1), eye(ord_a-1);-flip(a(2:end))];
    B = [zeros(ord_b-1,1), eye(ord_b-1);-flip(b(2:end))];

    G{3}  = [A, zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), B];

    Q{3}  = [o2w*diag(C_a), zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), o2v*diag(C_b)];

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

        %% FIltros de Kalman

        for m=1:NumModel
            x_{m} = G{m} * xa_{m};

            Pc{m} = (G{m} * Pa_{m} * G{m}') + ( Q{m} );

            z{m} = y(i,k) - (C{m}*x_{m});

            S{m} = C{m}*Pc{m}*C{m}' + R{m};

            K{m} = (Pc{m} * C{m}')/S{m};

            xo{m} = x_{m} + K{m} * (z{m});

            P{m} = (eye(ord_a+ord_b) - K{m} * C{m})* Pc{m};
        end

        %% Verossimilhança
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
        lkhd = zeros(NumModel,1);
        for m = 1:NumModel
            lkhd(m) = mvnpdf(z{m}, 0, S{m});
        end

        lkhd=lkhd+1e-50;
        nu = (Pr'*nu').*lkhd/( sum((Pr'*nu').*lkhd));
        nu = nu';
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
        %% Estimativa/Covariância Global

        %estimativa de x
        xt=sum( (nu.*cell2mat(xo))' )';

        %estimativa de P
        sm=[zeros(size(Pc{m}))];
        for m=1:NumModel
            sm=sm+nu(m) * (P{m} + (xo{m}-xt)*(xo{m}-xt)');
        end
        Pt=sm;
        clear sm

        xx(i,k)=xt(ord_a);
        out(i,k)=xt(ord_a-delay); %estimativa da saída com delay

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

        nur(i,k,1)=nu(1);
        nur(i,k,2)=nu(2);
        nur(i,k,3)=nu(3);
        j = j+1;

    end
    %...........................................................................................IMM
    nu=[mureinit(2:3)/(sum(mureinit(2:3))),0];
    xo=xreinit;
    P=Preinit;
end



ord=ord_a;

%% Adição de janelas sobrepostas 
OVERLAP_ADD
x=speech(1+delay:numel(t)+delay);

%Obtenção de nu para gráfico
out=nur(:,:,1);
nurr(:,1)=out(1:end);
out=nur(:,:,2);
nurr(:,2)=out(1:end);
out=nur(:,:,3);
nurr(:,3)=out(1:end);
end