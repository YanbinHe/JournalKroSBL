function [metrics] = svd_kroSBL_con2(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl)

% initialization
time = 0;

gamma1 = 1*ones(N,1); % prior on the gamma1
N1 = norm(gamma1);
gamma1 = gamma1/N1;
gamma2 = 1*ones(N,1); % prior on the gamma2
N2 = norm(gamma2);
gamma2 = gamma2/N2;
gamma3 = N1*N2*1*ones(N,1); % prior on the gamma3
gamma = kron(kron(gamma1,gamma2),gamma3);
x_re = 1*ones(N^3,1);% initialization

keep_list = 1:N^3;
keep_list1 = 1:N;
keep_list2 = 1:N;% absolute position
keep_list3 = 1:N;% absolute position

norms = 1;
itrs = 1;

r_absolute = 0;
r1 = 1e-3;
r2 = 1e-3;
r3 = 1e-3;
thres = 1e-4;

if func_ctrl.noisy_convergence
    R_max = 500;
    thres = 0;
end


while(itrs < R_max + 1 && norms > thres) % do the iteration
    itrs
    
    tic;
    
    gammal = zeros(N^3,1);
    gammal(keep_list) = gamma;
    gamma_old = gammal;
    
    if func_ctrl.noisy_convergence
        % Prune weights
        if min(abs(gamma1)) < r_absolute || min(abs(gamma2)) < r_absolute || min(abs(gamma3)) < r_absolute || min(abs(gamma1)) < r1*max(abs(gamma1)) || min(abs(gamma2)) < r2*max(abs(gamma2)) || min(abs(gamma3)) < r3*max(abs(gamma3)) 
            % if the first condition happens, it means that some block has to
            % be zero, if the second happens, it means that some entries in
            % each block should be zero.

            l1 = length(gamma1);
            l2 = length(gamma2);
            l3 = length(gamma3);   

            index1_dele = find(gamma1 < r1*max(abs(gamma1)));
            index1_sub = find(gamma1 < r_absolute);
            index1_dele = sort(unique([index1_dele;index1_sub]));

            index1 = 1:l1;
            index1(index1_dele) = [];% should keep relative position


            index2_dele = find(gamma2 < r2*max(abs(gamma2))); 
            index2_sub = find(gamma2 < r_absolute);
            index2_dele = sort(unique([index2_dele;index2_sub]));

            index2 = 1:l2;
            index2(index2_dele) = [];% should keep relative position

            index3_dele = find(gamma3 < r3*max(abs(gamma3))); 
            index3_sub = find(gamma3 < r_absolute);
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

            A = A(:,index);
            A1 = A1(:,index1);
            A2 = A2(:,index2);
            A3 = A3(:,index3);

            keep_list = keep_list(index);
            keep_list1 = keep_list1(index1);
            keep_list2 = keep_list2(index2);
            keep_list3 = keep_list3(index3);

            % prune gamma and corresponding entries in Sigma and mu
            gamma1 = gamma1(index1);
            gamma2 = gamma2(index2);
            gamma3 = gamma3(index3);

        end
    end
    
    % update the posterior mean and variance
    % compute the posterior
    [x_re,Sigma_x] = posterior_compute2(noise_var,gamma1,gamma2,gamma3,A1,A2,A3,A,y);

    Lambda = real(diag(Sigma_x)) + abs(x_re).^2;
    
    l1 = size(gamma1,1);
    l2 = size(gamma2,1);
    l3 = size(gamma3,1);
    
    mat1 = reshape(Lambda,l2*l3,l1);
    [mat1l,mat1v,mat1r] = svd(mat1');
    gamma1 = abs(mat1l(:,1));
    mat2 = reshape(abs(mat1v(1,1)*mat1r(:,1)),l3,l2);
    [mat2l,mat2v,mat2r] = svd(mat2');
    gamma2 = abs(mat2l(:,1));
    gamma3 = abs(mat2v(1,1)*mat2r(:,1));
    time = time + toc;   
    

        
    gammal = zeros(N^3,1);
    gamma = kron(gamma1,kron(gamma2,gamma3));
    gammal(keep_list) = gamma;
    
    norms = norm(gammal - gamma_old)/norm(gamma_old);
    norm(gammal - gamma_old)/norm(gamma_old)
    
    % compute error for convergence tracking
    x_rel = zeros(N^3,1);
    x_rel(keep_list) = x_re;
    
    errors(itrs) = norm(x_rel - x)/norm(x);
    
    itrs = itrs + 1;
end



if ~func_ctrl.noisy_convergence
    srr = recover_rate(x_rel,x);
    metrics={'error',errors;
             'vector',x_rel;
             'support recovery rate',srr;
             'time',time
             };
else
    metrics={'error',errors;
        'hyperparameter',[gamma1,gamma2,gamma3];
        };
end

end

