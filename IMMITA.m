function [x,murr]=IMMITA(n_ipt,fs,win_t,ord_a,ord_b,o2,noisesample,iter)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo


%% Inicialização do filtro
len = fix(win_t * fs);


%% Cortando a entrada nas janelas que vamos trabalhar
len = fix(win_t * fs);
len_ovrlp=round(len/2);

y=buffer(n_ipt,len,len_ovrlp);
n=size(y);
delay=ord_a-1;


%% Inicialização do filtro

NumModel=3;
H_ = [zeros(1,ord_a-1),1];
H{1} = [zeros(1,ord_a-1),1,zeros(1,ord_b-1),1];
R{1} = o2;
P{1} = R{1} * eye(ord_a+ord_b);

x = zeros(1,length(n_ipt));
xa = zeros(ord_a+ord_b,1);

[b, o2v] = lpc(noisesample.*hamming(numel(noisesample)),ord_b);

j = 1;

MarkovModel
mu=[0.75,0.25,0.0];
xo{1}=xa;

for ii=2:NumModel
    P{ii}=P{1};
    xo{ii}=xo{1};
    H{ii}=H{1};
    R{ii}=R{1};
end

xt=xo{1};
Pt=P{1};

[A{3},o2w{3}] = KFA_LPC(y(:,1),o2,ord_a,ord_b,iter-1,b,o2v,xt,Pt); %Obtenção do modelo através de um filtro de Kalman Aumentado
A{2}=A{3};o2w{2}=o2w{3};
%% Processamento do sinal
for k = 1:n(2)

    A{1}=A{2};o2w{1}=o2w{2};
    A{2}=A{3};o2w{2}=o2w{3};
    [A{3},o2w{3}] = KFA_LPC(y(:,k+1*(k<n(2))),o2,ord_a,ord_b,iter-1,b,o2v,xt,Pt); %Obtenção do modelo através de um filtro de Kalman Aumentado

    %...........................................................................................IMM
    for i = 1:len
        %% Reinicialização

        mu_=(Pr.*mu./(mu*Pr)')';

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
            x_{m} = A{m} * xa_{m};

            Pc{m} = (A{m} * Pa_{m} * A{m}') + ( o2w{m} );

            z{m} = y(i,k) - (H{m}*x_{m});

            S{m} = H{m}*Pc{m}*H{m}' + R{m};

            K{m} = (Pc{m} * H{m}')/S{m};

            xo{m} = x_{m} + K{m} * (z{m});

            P{m} = (eye(ord_a+ord_b) - K{m} * H{m})* Pc{m};
        end

        %% Verossimilhança


        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
        lkhd = zeros(NumModel,1);
        for m = 1:NumModel
            lkhd(m) = mvnpdf(z{m}, 0, S{m});
        end

        lkhd=lkhd+1e-50;
        mu = (Pr'*mu').*lkhd/( sum((Pr'*mu').*lkhd));
        mu = mu';
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
        %% Estimativa/Covariância Global

        %estimativa de x
        xt=sum( (mu.*cell2mat(xo))' )';

        %estimativa de P
        sm=[zeros(size(Pc{m}))];
        for m=1:NumModel
            sm=sm+mu(m) * (P{m} + (xo{m}-xt)*(xo{m}-xt)');
        end
        Pt=sm;
        clear sm

        xx(i,k)=xt(ord_a);
        out(i,k)=xt(ord_a-delay);

        %   Armazenamento de variáveis para transição para janela futura
        if i==len-len_ovrlp
            xreinit{1} = xo{2};
            xreinit{2} = xo{3};
            xreinit{3} = xt;
            Preinit{1} = P{2};
            Preinit{2} = P{3};
            Preinit{3} = Pt;
            mureinit = mu;
        end

        mur(i,k,1)=mu(1);
        mur(i,k,2)=mu(2);
        mur(i,k,3)=mu(3);
        j = j+1;

    end
    %...........................................................................................IMM
    mu=[mureinit(2:3),mureinit(1)];
    xo=xreinit;
    P=Preinit;
end


%% Adição de janelas sobrepostas
ord=ord_a;
OVERLAP_ADD
x=speech(1+delay:numel(t)+delay);

out=mur(:,:,1);
murr(:,1)=out(1:end);
out=mur(:,:,2);
murr(:,2)=out(1:end);
out=mur(:,:,3);
murr(:,3)=out(1:end);
end