function [metrics] = svd_kroSBL_con_nsq(y,A1,A2,A3,A,N,R_max,x,func_ctrl)

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
num_M = length(y);

keep_list = 1:N^3;
keep_list1 = 1:N;
keep_list2 = 1:N;% absolute position
keep_list3 = 1:N;% absolute position
noise_var = 1e-1*var(y);

norms = 1;
itrs = 1;

r_absolute = 0;
r1 = 1e-3;
r2 = 1e-3;
r3 = 1e-3;
thres = 1e-4;


R_max = 200;
thres = 0;



while(itrs < R_max + 1 && norms > thres) % do the iteration
    itrs

    tic;
    
    gammal = zeros(N^3,1);
    gammal(keep_list) = gamma;
    gamma_old = gammal;
    

    gammad = kron(gamma1,kron(gamma2,gamma3));
    
    % update the posterior mean and variance
    % compute the posterior
    [x_re,Sigma_x] = posterior_compute2(noise_var,gamma1,gamma2,gamma3,A1,A2,A3,A,y);

    Lambda = real(diag(Sigma_x)) + abs(x_re).^2;Sigma_x_diag = real(diag(Sigma_x));
    
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

    % update noise variance
    temp = sum(1 - Sigma_x_diag./gammad);
    noise_var = (norm(y - A*x_re)^2 + noise_var*temp)/num_M;
    
    time = time + toc;   
    

        
    gammal = zeros(N^3,1);
    gamma = kron(gamma1,kron(gamma2,gamma3));
    gammal(keep_list) = gamma;
    
    norms = norm(gammal - gamma_old)/norm(gamma_old);
    norm(gammal - gamma_old)/norm(gamma_old);
    
    % compute error for convergence tracking
    
    errors(itrs) = (norm(x_re - x)/norm(x))^2;(norm(x_re - x)/norm(x))^2
    
    itrs = itrs + 1;
end

metrics={'error',errors;
    'hyperparameter',[gamma1,gamma2,gamma3];
    };

end

