function [a]=CorrV(x,n)
    C=xcorr(x,n,'biased');
    a=C(n+1:end);
end
