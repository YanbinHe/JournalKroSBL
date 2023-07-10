function [rankl] = svd_kroSBL_noiseless(y1,y2,y3,A1,A2,A3,N,R_max,x)

% initialization

gamma1 = 1*ones(N,1); % prior on the gamma1
N1 = norm(gamma1);
gamma1 = gamma1/N1;
gamma2 = 1*ones(N,1); % prior on the gamma2
N2 = norm(gamma2);
gamma2 = gamma2/N2;
gamma3 = N1*N2*1*ones(N,1); % prior on the gamma3
gamma = kron(gamma1,kron(gamma2,gamma3));

x_re = ones(N^3,1);

norms = 1;
itrs = 1;

thres = 0;

while(itrs < R_max + 1 && norms > thres) % do the iteration
    itrs
    
    gamma_old = gamma;
    
    % update the posterior mean and variance
    % compute the posterior
    [x_re,Sigma_x] = posterior_compute_noiseless(gamma1,gamma2,gamma3,A1,A2,A3,y1,y2,y3);

    Lambda = real(diag(Sigma_x)) + abs(x_re).^2;
    
    [rankl(1,itrs),rankl(2,itrs)] = rank_compute(Lambda,N);

    mat1 = reshape(Lambda,N^2,N);
    [mat1l,mat1v,mat1r] = svd(mat1');
    gamma1 = abs(mat1l(:,1));
    mat2 = reshape(abs(mat1v(1,1)*mat1r(:,1)),N,N);
    [mat2l,mat2v,mat2r] = svd(mat2');
    gamma2 = abs(mat2l(:,1));
    gamma3 = abs(mat2v(1,1)*mat2r(:,1));
  
    % re-estimate
    gamma = kron(gamma1,kron(gamma2,gamma3));
%     [x_re,~] = posterior_compute_noiseless(gamma1,gamma2,gamma3,A1,A2,A3,y1,y2,y3);
        
    
    norms = norm(gamma - gamma_old)/norm(gamma_old) 
    
%     errors = norm(x_re - x)/norm(x)

    itrs = itrs + 1;
end
end

