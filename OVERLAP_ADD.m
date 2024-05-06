
if len_ovrlp>0
wgh=hann(len);
wgth=wgh(1:len_ovrlp);
siz=size(out);

speech=[out(len_ovrlp+1:end-len_ovrlp,1)];

for k=2:siz(2)
    speech = [speech; out(1:len_ovrlp,k).*wgth + out(end-len_ovrlp+1:end,k-1).*(1-wgth)];
    speech = [speech; out(len_ovrlp+1:end-len_ovrlp,k)];
end
speech=[speech;out(end-len_ovrlp+1:end,end);zeros(ord,1)];
else
    speech=[out(1:end)';zeros(ord,1)];
end

