function [x,out]=KFIT(n_ipt,fs,win_t,ord,o2w,iter)

t = (0:1/fs:(length(n_ipt)-1)/fs)';              %Definição de tempo


%% Inicialização do filtro
len = fix(win_t * fs);
len_ovrlp=round(len/2);

% Cortando a entrada nas janelas que vamos trabalhar
n=floor(length(n_ipt)/len);
y=buffer(n_ipt,len,len_ovrlp);

H = [zeros(1,ord-1),1];
R = o2w;
wndw=hamming(len);

P = R * eye(ord);
x = zeros(1,length(n_ipt));
xo = zeros(ord,1);
xo(end)= n_ipt(1);
larg=size(y);
delay=ord-1;

%% Processamento do sinal
for k=1:larg(2)

    % Obtenção do modelo para primeira iteração
    [lpccoef,o2w] = lpc(y(:,k).*wndw, ord);

    %   Armazenamento de variáveis para reinicio de janelas
    xreinit1 = xo;
    Preinit1 = P;


    for l = 1:iter
        lpcrec=lpccoef;
        A = [zeros(ord-1,1), eye(ord-1); -flip(lpccoef(2:end))];
        for i = 1:len
            %   Filtro de Kalman
            x_ = A * xo;
            Pc = (A * P * A') + (o2w * diag(H));
            K = (Pc * H')/((H * Pc * H') + R);
            xo = x_ + K * (y(i,k) - (H*x_));
            P = (eye(ord) - K * H) * Pc;

            %   Estimativa da saída com delay
            xx(i,k)=xo(end);
            out(i,k)=xo(end-delay);

            % obtenção de gráfico de K_q,k
            Kr(k,i)=K(end);

            %   Armazenamento dde x_k' e P_k'
            if i==len-len_ovrlp
                xreinit2 = xo;
                Preinit2 = P;
            end

        end
        % Obtenção do modelo para iteração l+1
        [lpccoef,o2w] = lpc(xx(:,k).*wndw, ord);
        xo = xreinit1*(l<iter) + xo*(l==iter);
        P = Preinit1*(l<iter) +P*(l==iter);
    end
    xo=xreinit2;
    P=Preinit2;
end

%% Adição de janelas sobrepostas 
OVERLAP_ADD
x=speech(1+delay:numel(t)+delay);

% Adição de janelas sobrepostas (K_p,k)
% out=Kr;
% OVERLAP_ADD
% Kr=speech(1+delay:numel(t)+delay);
end