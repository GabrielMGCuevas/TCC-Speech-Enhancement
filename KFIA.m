function [x,cnoise]=KFITA(n_ipt,fs,win_t,ord_a,ord_b,o2,noisesample,iter)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo
x = zeros(1,length(n_ipt));
cnoise = zeros(1,length(n_ipt));


%% Inicialização do filtro

len = fix(win_t * fs);
len_ovrlp=round(len/2);
% Cortando a entrada nas janelas que vamos trabalhar
n=floor(length(n_ipt)/len);
y=buffer(n_ipt,len,len_ovrlp);

wndw=hamming(len);
wndw2=hamming(numel(noisesample));

R = o2;
P = R * eye(ord_b);
x = zeros(1,length(n_ipt));
xo = zeros(ord_b,1);
larg=size(y);
delay=ord_a-1;


%% Modelo do ruído

[b, o2u] =lpc(noisesample.*hamming(numel(noisesample)),ord_b);

C= [zeros(1,ord_a-1),1,zeros(1,ord_b-1),1];
C_=[zeros(1,ord_a-1),1];
xo = zeros(ord_a+ord_b,1);
R = o2;
B = [zeros(ord_b-1,1), eye(ord_b-1);...
    -flip(b(2:end))];

o2u = o2u*diag([zeros(1,ord_b-1),1]);

P = R * eye(ord_a+ord_b);

%% Processing.
for k=1:larg(2)

    % Obtenção do modelo da primeira iteração através da aplicação do WF
    [a,o2w] = lpc(filter(b,1,y(:,k)).*wndw, ord_a);  
    o2w=sum(CorrV(y(:,k).*wndw,ord_a).*a')    -   o2u(end,end);
    xreinit1 = xo;
    Preinit1 = P;

    for l = 1:iter


        A  = [zeros(ord_a-1,1), eye(ord_a-1);...
            -flip(a(2:end))];

        G  = [A, zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), B];

        Q  = [o2w*diag(C_), zeros(ord_a,ord_b);...
            zeros(ord_b,ord_a), o2u];

        for i = 1:len
            % Filtro de Kalman
            x_ = G * xo;
            Pc = (G * P * G') +  Q;
            K = (Pc * C')/((C * Pc * C') + R);
            xo = x_ + K * (y(i,k) - (C*x_));
            P = (eye(length(G)) - K * C) * Pc;


            cnoise(j)=xo(end);

            xx(i,k)=xo(ord_a);
            out(i,k)=xo(ord_a-delay); %estimativa da saída com delay

            if i==len_ovrlp
                xreinit2 = xo;
                Preinit2 = P;
            end
        end

        [a,o2w] = lpc(xx(:,k).*wndw, ord_a); % Obtenção do modelo da n-ésima iteração

        j = j-(len)*(l<iter);
        xo = xreinit1*(l<iter) + xo*(l==iter);
        P = Preinit1*(l<iter) + P*(l==iter);

    end

    j=j-(len_ovrlp+1);
    xo=xreinit2;
    P=Preinit2;
end

ord=ord_a;
%% Adição de janelas sobrepostas 
OVERLAP_ADD
x=speech(1+delay:numel(t)+delay);


end