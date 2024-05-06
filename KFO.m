function [s,xx,Kr]=KFO(n_ipt,fs,win_t,ord,o2,iter,dataC)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo


%% Inicialização do filtro
len = fix(win_t * fs);
len_ovrlp=round(len/2);

% Cortando a entrada nas janelas que vamos trabalhar
n=floor(length(n_ipt)/len);
y=buffer(n_ipt,len,len_ovrlp);
yc=buffer(dataC,len,len_ovrlp);

% Definição de parêmtros do filtro
C = [zeros(1,ord-1),1];
R = o2;
wndw=hamming(len);

P = R * eye(ord);
s = zeros(1,length(n_ipt));
xo = zeros(ord,1);
xo(end)= n_ipt(1);
larg=size(y);

delay=ord-1;
%% Processamento do sinal
for k=1:larg(2)

    %   Obtenção do modelo
    [lpccoef,o2w] = lpc(yc(:,k).*wndw, ord);

    for l = 1:1
        lpcrec=lpccoef;
        A = [zeros(ord-1,1), eye(ord-1); -flip(lpccoef(2:end))];

        for i = 1:len
             %   Filtro de Kalman
            x_ = A * xo;
            Pc = (A * P * A') + (o2w * diag(C));
            K = (Pc * C')/((C * Pc * C') + R);
            xo = x_ + K * (y(i,k) - (C*x_));
            P = (eye(ord) - K * C) * Pc;

            %   Estimativa da saída com delay
            out(i,k)=xo(end-delay);

            % obtenção de gráfico de K_q,k
            Kr(i,k)=K(end);

             %   Armazenamento dde x_k' e P_k'
            if i==len-len_ovrlp
                xreinit = xo;
                Preinit = P;
            end
        end
    end
    xo=xreinit;
    P=Preinit;
end

%% Adição de janelas sobrepostas 
OVERLAP_ADD
s=speech(1+delay:numel(t)+delay);

% Adição de janelas sobrepostas (K_p,k)
out=Kr;
OVERLAP_ADD
Kr=speech(1+delay:numel(t)+delay);

end