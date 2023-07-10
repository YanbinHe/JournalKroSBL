function [metrics] = OLS(y,A,N,support,x)
% oracle least square

xset = pinv(A(:,support))*y;
x_re = zeros(N^3,1);
x_re(support) = xset;
errorols = norm(x_re - x)/norm(x) 

metrics={'error',errorols;
         'vector',x_re;
         'support recovery rate',1;
         'time',0
         };

end

