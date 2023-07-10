function [errorl,time,mu_x,gammare,H_re] = kroSBL_pcsvd(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,irs_pattern,A_irs_a,A_irs_d,H1,H2,SNR,thres,overheadi,alpha,beta)

time = 0;

A2_ori = A2; % a copy of the dictionary

alpha_diag = alpha*ones(Res1-1,1);
beta_diag = beta*ones(Res1-1,1);
Calpha = eye(Res1) + diag(alpha_diag,-1) + diag(alpha_diag,1);
Cbeta = eye(Res1) + diag(beta_diag,-1) + diag(beta_diag,1);

gamma1 = 1*ones(Res1,1); % prior on the gamma1
gamma1t = Calpha*gamma1;
gamma2 = 1*ones(Res1,1); % prior on the gamma2
gamma3 = 1*ones(Res1,1); % prior on the gamma3
gamma3t = Cbeta*gamma3;

mu_x = 1*ones(Res1^3,1);% initialization

norm1 = 1;
itr1 = 1;

while(itr1 < numItr + 1 && norm1 > thres) % do the iteration
    itr1
    
    tic;
    
    mu_x_old = mu_x;
     
    % update the posterior mean and variance
    % compute the posterior
    [mu_x,Sigma_x] = posterior_compute2(noise_var,gamma1t,gamma2,gamma3t,H_p1,H_p2,A2,H,y);

    Lambda = real(diag(Sigma_x)) + abs(mu_x).^2;
    
    l1 = size(gamma1t,1);
    l2 = size(gamma2,1);
    l3 = size(gamma3t,1);
    
    mat1 = reshape(Lambda,l2*l3,l1);
    [mat1l,mat1v,mat1r] = svd(mat1');
    gamma1t = abs(mat1l(:,1));
    gamma1 = lsqnonneg(Calpha,gamma1t);
    gamma1t = Calpha*gamma1;
    
    mat2 = reshape(abs(mat1v(1,1)*mat1r(:,1)),l3,l2);
    [mat2l,mat2v,mat2r] = svd(mat2');
    gamma2 = abs(mat2l(:,1));
    
    gamma3t = abs(mat2v(1,1)*mat2r(:,1));
    gamma3 = lsqnonneg(Cbeta,gamma3t);
    gamma3t = Cbeta*gamma3;
    time = time + toc;   
    
    % re-estimate
    [mu_x,~] = posterior_compute2(noise_var,gamma1t,gamma2,gamma3t,H_p1,H_p2,A2,H,y);
        
    gammare = kron(gamma1t,kron(gamma2,gamma3t));
    
    norm1 = norm(mu_x - mu_x_old)/norm(mu_x_old)
    
    % computing the error for each IRS pattern and averaging them all  
    for i = 1:K2
        Htt = irs_pattern(:,i).'*kr(A_irs_a.',A_irs_d').';
        Ht(:,:,i) = Htt(:,1:Res1);
        H_re(:,:,i) = kron(kron(Ht(:,:,i),conj(A1)),A2_ori)*mu_x;
        Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
        nmse(itr1,i) = norm(H_re(:,:,i) - Htrue(:,:,i),'fro')/norm(Htrue(:,:,i),'fro');
    end
    
    error = sum(nmse(itr1,:))/K2
    errorl(itr1) = error;
    
    itr1 = itr1 + 1;
end
errorl = errorl(end);
end

