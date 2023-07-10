function [rankl] = am_kroSBL_noiseless(y1,y2,y3,A1,A2,A3,N,R_max,x)
gamma1 = 1*ones(N,1); % prior on the gamma1
N1 = norm(gamma1);
gamma1 = gamma1/N1;
gamma2 = 1*ones(N,1); % prior on the gamma2
N2 = norm(gamma2);
gamma2 = gamma2/N2;
gamma3 = N1*N2*1*ones(N,1); % prior on the gamma3
gamma = kron(gamma1,kron(gamma2,gamma3));

x_re = ones(N^3,1);

norma = 1;
itra = 1;

thres_inner = 1e-10;
thres = 0;

while(itra < R_max + 1 && norma > thres) % do the iteration
    itra
    
    gamma_old = gamma;
    % update the posterior mean and variance
    % compute the posterior
    [x_re,Sigma_x] = posterior_compute_noiseless(gamma1,gamma2,gamma3,A1,A2,A3,y1,y2,y3);

    Lambda = real(diag(Sigma_x)) + abs(x_re).^2;   
    
    [rankl(1,itra),rankl(2,itra)] = rank_compute(Lambda,N);
    
    normi1 = 1;

    while(normi1 > thres_inner)
        
        gamma_old_inner = kron(kron(gamma1,gamma2),gamma3);

        % update 1 and projection
        gamma1 = N^-2*(kron(eye(N),kron(gamma2.^-1,gamma3.^-1)))'*Lambda;
        % update 2 and projection
        gamma2 = N^-2*(kron(gamma1.^-1,kron(eye(N),gamma3.^-1)))'*Lambda;
        % update 3
        gamma3 = N^-2*(kron(gamma1.^-1,kron(gamma2.^-1,eye(N))))'*Lambda;
        n1 = norm(gamma1);
        n2 = norm(gamma2);
        gamma1 = gamma1/n1; 
        gamma2 = gamma2/n2;
        gamma3 = gamma3*n1*n2;
        
        normi1 = norm(kron(kron(gamma1,gamma2),gamma3) - gamma_old_inner);        
    end
 
    gamma = kron(gamma1,kron(gamma2,gamma3));
    % re-estimate
%     [x_re,~] = posterior_compute_noiseless(gamma1,gamma2,gamma3,A1,A2,A3,y1,y2,y3);
    
    norma = norm(gamma - gamma_old)/norm(gamma_old) 
    
%     errora = norm(x_re - x)/norm(x) 
    
    itra = itra + 1;
end
end

