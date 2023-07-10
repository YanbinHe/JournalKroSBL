function s = pilot_gen(d1,d2)
% construct a QPSK source matrix
s = zeros(d1*d2,1);
for i = 1:d1*d2
    seed = rand(1);
    if seed>=0&&seed<0.25
        s(i)=(1+1i)/sqrt(2);
    elseif seed>=0.25&&seed<0.5
        s(i)=(1-1i)/sqrt(2);
    elseif seed>=0.5&&seed<0.75
        s(i)=(-1+1i)/sqrt(2);
    elseif seed>=0.75&&seed<=1
        s(i)=(-1-1i)/sqrt(2);
    end
end

s = reshape(s,d1,d2);
end