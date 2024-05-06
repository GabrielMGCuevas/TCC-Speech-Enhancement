function [a,Q,xo,P] = KF_LPC2(y,o2,ord,iter,xo,P)

        len=numel(y);
        wndw=hamming(len);
        

        H = [zeros(1,ord-1),1];  
        R = o2;

        [a,Q] = lpc(y.*wndw, ord);

        for oo = 1:iter
                A = [zeros(ord-1,1), eye(ord-1); -flip(a(2:end))];
        
                        for i = 1:len
                            x_ = A * xo;
                            Pc = (A * P * A') + (Q * diag(H));
                            K = (Pc * H')/((H * Pc * H') + R);
                            xo = x_ + K * (y(i) - (H*x_));
                            P = (eye(ord) - K * H) * Pc;
                            xx(i)=xo(end);
                        end     
                [a,Q] = lpc(xx.*wndw', ord);
        end

end

