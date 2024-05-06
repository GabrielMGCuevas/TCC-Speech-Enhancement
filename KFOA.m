function [x]=KFOA(n_ipt,fs,win_t,ord_a,ord_b,o2,noisesample,iter,c_ipt)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo
x = zeros(1,length(n_ipt));
cnoise = zeros(1,length(n_ipt));

%% Inicialização do filtro

len = fix(win_t * fs);
len_ovrlp=round(len/2);
% Cortando a entrada nas janelas que vamos trabalhar
n=floor(length(n_ipt)/len);

y=buffer(n_ipt,len,len_ovrlp);
yc=buffer(c_ipt,len,len_ovrlp);

wndw=hamming(len);
wndw2=hamming(numel(noisesample));

% Definição de parêmtros do filtro
R = o2;
C = R * eye(ord_b);
x = zeros(1,length(n_ipt));
xo = zeros(ord_b,1);
j = 1; % para criação de saída
larg=size(y);
delay=ord_a-1;



%% Modelo do ruído

H= [zeros(1,ord_a-1),1,zeros(1,ord_b-1),1];
H_a=[zeros(1,ord_a-1),1];
H_b=[zeros(1,ord_b-1),1];
xo = zeros(ord_a+ord_b,1);
R = o2;
C = R * eye(ord_a+ord_b);

%% Processamento do sinal
for k=1:larg(2)

   %   Obtenção do modelo
    [a,o2w] = lpc(yc(:,k).*wndw, ord_a);
    [b,o2v] = lpc((y(:,k)-yc(:,k)).*wndw, ord_b);

    for l = 1:1

            %   Filtro de Kalman
        A  = [zeros(ord_a-1,1), eye(ord_a-1);-flip(a(2:end))];
        B = [zeros(ord_b-1,1), eye(ord_b-1);-flip(b(2:end))];
        G  = [A, zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), B];
        Q  = [o2w*diag(H_a), zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), o2v*diag(H_b)];

        for i = 1:len

            x_ = G * xo;
            Pc = (G * C * G') +  Q;
            K = (Pc * H')/((H * Pc * H') + R);
            xo = x_ + K * (y(i,k) - (H*x_));
            C = (eye(length(G)) - K * H) * Pc;


            %   Estimativa da saída com delay
            xx(i,k)=xo(ord_a);
            out(i,k)=xo(ord_a-delay);


             %   Armazenamento dde x_k' e P_k'
            if i==len_ovrlp
                xreinit2 = xo;
                Preinit2 = C;
            end
        end
    end

    xo=xreinit2;
    C=Preinit2;
end

ord=ord_a;
%% Adição de janelas sobrepostas 
OVERLAP_ADD
x=speech(1+delay:numel(t)+delay);


end