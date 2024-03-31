function [rankl] = classicSBL_noiseless(y,A,N,R_max,x)

time = 0;

gammac = 1*ones(N^3,1); % initialize the prior
x_re = ones(N^3,1);

normc = 1;
itrc = 1;

thres = 0;

while(itrc < R_max + 1 && normc > thres) % do the iteration
    itrc
    
    gamma_old = gammac;
    % update the posterior
    Gammac = diag(gammac);
    Gammacsqrt = diag(gammac.^0.5);
    
    temp = Gammacsqrt*pinv(A*Gammacsqrt);
    Sigma_x = Gammac - temp*A*Gammac;
    x_re = temp * y;
    
    % updating gamma
    Sigma_x_diag = real(diag(Sigma_x));
    gammac = Sigma_x_diag + abs(x_re).^2;
    
    [rankl(1,itrc),rankl(2,itrc)] = rank_compute(gammac,N);
    
    % re-estimate
%     Gammacsqrt = diag(gammac.^0.5);
%     
%     temp = Gammacsqrt*pinv(A*Gammacsqrt);
%     x_re = temp * y;

    normc = norm(gammac-gamma_old)/norm(gamma_old)
    
%     errorc = norm(x_re - x)/norm(x)

    itrc = itrc + 1; 
end
end

