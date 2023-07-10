function s = qpsk_sym(sym)
% construct a QPSK source sequence s
s = zeros(sym,1);
for i = 1:sym
    seed = rand(1);
    if seed>=0&&seed<0.25
        s(i)=(1+1i);
    elseif seed>=0.25&&seed<0.5
        s(i)=(1-1i);
    elseif seed>=0.5&&seed<0.75
        s(i)=(-1+1i);
    elseif seed>=0.75&&seed<=1
        s(i)=(-1-1i);
    end
end
end