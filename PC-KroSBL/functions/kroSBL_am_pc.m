function [errorl,time,mu_x,gammare,H_re] = kroSBL_am_pc(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,irs_pattern,A_irs_a,A_irs_d,H1,H2,SNR,thres,overheadi,alpha,beta)
% this function is formulated using variance
time = 0;

A2_ori = A2; % a copy of the dictionary
global thres_inner;

alpha_diag = alpha*ones(Res1-1,1);
beta_diag = beta*ones(Res1-1,1);
Calpha = eye(Res1) + diag(alpha_diag,-1) + diag(alpha_diag,1);
Cbeta = eye(Res1) + diag(beta_diag,-1) + diag(beta_diag,1);

mu_x_est = 1*ones(Res1^3,1);% initialization

gamma1 = 1*ones(Res1,1); % prior on the gamma1
gamma1t = Calpha*gamma1;
gamma2 = 1*ones(Res1,1); % prior on the gamma2
gamma3 = 1*ones(Res1,1); % prior on the gamma3
gamma3t = Cbeta*gamma3;
gamma = kron(gamma1t,kron(gamma2,gamma3t));

keep_list = 1:Res1^3;
keep_list1 = 1:Res1;
keep_list2 = 1:Res1;% absolute position
keep_list3 = 1:Res1;% absolute position

norm1 = 1;
itr1v = 1;

r_absolute = 1e-4;
r1 = 1e-2;
r2 = 1e-1;
r3 = 1e-2;


while(itr1v < numItr + 1 && norm1 > thres) % do the iteration
    itr1v
    
    tic;
    
    gammare = zeros(Res1^3,1);
    gammare(keep_list) = gamma;
    gamma_old = gammare;
 
    % Prune weights
    if min(abs(gamma1t)) < r_absolute || min(abs(gamma2)) < r_absolute || min(abs(gamma3t)) < r_absolute || min(abs(gamma1t)) < r1*max(abs(gamma1t)) || min(abs(gamma2)) < r2*max(abs(gamma2)) || min(abs(gamma3t)) < r3*max(abs(gamma3t)) 
        % if the first condition happens, it means that some block has to
        % be zero, if the second happens, it means that some entries in
        % each block should be zero.
        
        l1 = length(gamma1t);
        l2 = length(gamma2);
        l3 = length(gamma3t);   
        
        index1_dele = find(gamma1t < r1*max(abs(gamma1t)));
        index1_sub = find(gamma1t < r_absolute);
        index1_dele = sort(unique([index1_dele;index1_sub]));
        
        index1 = 1:l1;
        index1(index1_dele) = [];% should keep relative position
        
        
        index2_dele = find(gamma2 < r2*max(abs(gamma2))); 
        index2_sub = find(gamma2 < r_absolute);
        index2_dele = sort(unique([index2_dele;index2_sub]));
        
        index2 = 1:l2;
        index2(index2_dele) = [];% should keep relative position
        
        index3_dele = find(gamma3t < r3*max(abs(gamma3t))); 
        index3_sub = find(gamma3t < r_absolute);
        index3_dele = sort(unique([index3_dele;index3_sub]));
        
        index3 = 1:l3;
        index3(index3_dele) = [];% should keep relative position
        
        index = [];%zeros(l1*l2,1);
        
        for i =1:length(index1)
            for j = 1:length(index2)
                for k = 1:length(index3)
                    index = [index;(((index1(i)-1)*l2+index2(j)-1)*l3+index3(k))];
                end
            end
        end
        
        H = H(:,index);
        H_p1 = H_p1(:,index1);
        H_p2 = H_p2(:,index2);
        A2 = A2(:,index3);
        
        keep_list = keep_list(index);
        keep_list1 = keep_list1(index1);
        keep_list2 = keep_list2(index2);
        keep_list3 = keep_list3(index3);

        % prune gamma and corresponding entries in Sigma and mu
        gamma1t = gamma1t(index1);
        gamma2 = gamma2(index2);
        gamma3t = gamma3t(index3);
        
        Calpha = Calpha(index1,index1);
        Cbeta = Cbeta(index3,index3);

    end
    
    % update the posterior mean and variance
    % compute the posterior
    [mu_x,Sigma_x] = posterior_compute2(noise_var,gamma1t,gamma2,gamma3t,H_p1,H_p2,A2,H,y);
    
    Lambda = real(diag(Sigma_x)) + abs(mu_x).^2;   
        
    l1 = size(gamma1t,1);
    l2 = size(gamma2,1);
    l3 = size(gamma3t,1);

    normi1 = 1;
    itr_am = 1;

    while(normi1 > thres_inner && itr_am < 10)
        
        gamma_old_inner = kron(gamma1t,kron(gamma2,gamma3t));
        % update 1
        gamma1 = (lsqnonneg(Calpha,(l2*l3)^-1*((kron(eye(l1),kron(gamma2.^-1,gamma3t.^-1)))'*Lambda)));
        gamma1t = 1e-5 + Calpha*gamma1;
        
        gamma2 = (l1*l3)^-1*((kron(gamma1t.^-1,kron(eye(l2),gamma3t.^-1)))'*Lambda);
         
        gamma3 = (lsqnonneg(Cbeta,(l1*l2)^-1*((kron(gamma1t.^-1,kron(gamma2.^-1,eye(l3))))'*Lambda)));
        gamma3t = 1e-5 + Cbeta*gamma3;
        
        n1 = norm(gamma1t);
        gamma1t = gamma1t/n1; 
        n2 = norm(gamma2);
        gamma2 = gamma2/n2;
        gamma3t = gamma3t*n1*n2;

        normi1 = norm(kron(kron(gamma1t,gamma2),gamma3t) - gamma_old_inner)/norm(gamma_old_inner);
        itr_am = itr_am + 1;
    end
    
    time = time + toc;   
    gamma = kron(kron(gamma1t,gamma2),gamma3t);
    gammare = zeros(Res1^3,1);
    gammare(keep_list) = gamma;
    
    norm1 = norm(gammare - gamma_old)/norm(gamma_old);
    norm(gammare - gamma_old)/norm(gamma_old)
    
    itr1v = itr1v + 1;
end
% re-estimate
[mu_x,~] = posterior_compute2(noise_var,gamma1t,gamma2,gamma3t,H_p1,H_p2,A2,H,y);

mu_x_est = zeros(Res1^3,1);
mu_x_est(keep_list) = mu_x;

% computing the error for each IRS pattern and averaging them all
for i = 1:K2
    Htt = irs_pattern(:,i).'*kr(A_irs_a.',A_irs_d').';
    Ht(:,:,i) = Htt(:,1:Res1);
    H_re(:,:,i) = kron(kron(Ht(:,:,i),conj(A1)),A2_ori)*mu_x_est;
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    nmse(i) = norm(H_re(:,:,i) - Htrue(:,:,i),'fro')/norm(Htrue(:,:,i),'fro');
end

error = sum(nmse(:))/K2
errorl = error;
end

