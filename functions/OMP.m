function [metrics] = OMP(noise_var,y,A,N,x)

time = 0;

tic;

residual = y;
set = [];
x_re = zeros(N^3,1);

while norm(residual) > 1*sqrt(noise_var)*sqrt(length(y))

    inner_pro = abs(A'*residual);
    [~,in] = max(inner_pro);
    set = [set in];
    
    Aset = A(:,set);
    xset = pinv(Aset)*y;
    
    x_re(set) = xset;
    residual = y-A*x_re;
end

time = time + toc;
erroro = norm(x_re - x)/norm(x)

srr = recover_rate(x_re,x);
metrics={'error',erroro;
         'vector',x_re;
         'support recovery rate',srr;
         'time',time
         };
end

