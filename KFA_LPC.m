function [G,Q,xx] = KFA_LPC(y,o2,ord,ord2,iter,b,o2u,xo,P)

        len=numel(y);
        wndw=hamming(len);
        
        B = [zeros(ord2-1,1), eye(ord2-1);...
            -flip(b(2:end))];

        C_ = [zeros(1,ord-1),1];  
        H= [zeros(1,ord-1),1,zeros(1,ord2-1),1];
        Q_=o2u*diag([zeros(1,ord2-1),1]);

        R = o2;    
        x = zeros(1,len); 

    [a,o2w] = lpc(filter(b,1,y(:)).*wndw, ord);
    o2w=sum(CorrV(y.*wndw,ord).*a')    -   Q_(end,end); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        for oo = 1:iter
               A  = [zeros(ord-1,1), eye(ord-1);...
                      -flip(a(2:end))];
            
                G  = [A, zeros(ord,ord2);...
                    zeros(ord2,ord), B];
            
                Q  = [o2w*diag(C_), zeros(ord,ord2);...
                    zeros(ord2,ord), Q_];   
        
                        for i = 1:len
                            x_ = G * xo;
                            Pc = (G * P * G') +  Q;
                            K = (Pc * H')/((H * Pc * H') + R);
                            xo = x_ + K * (y(i) - (H*x_));
                            P = (eye(ord+ ord2) - K * H) * Pc;

                            xx(i)=xo(ord);
                            rr(i)=xo(end);
                        end     
        
                [a,o2w] = lpc(xx.*wndw', ord);
                a=a;
        end

                A  = [zeros(ord-1,1), eye(ord-1);...
                      -flip(a(2:end))];
            
                G  = [A, zeros(ord,ord2);...
                    zeros(ord2,ord), B];
        
                Q  = [o2w*diag(C_), zeros(ord,ord2);...
                    zeros(ord2,ord), Q_];

end

