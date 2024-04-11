function [metrics] = am_kroSBL_con3(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl)

metrics = {};
time = 0;
gamma1 = 1*ones(N,1); % prior on the gamma1
N1 = norm(gamma1);
gamma1 = gamma1/N1;
gamma2 = 1*ones(N,1); % prior on the gamma2
N2 = norm(gamma2);
gamma2 = gamma2/N2;
gamma3 = N1*N2*1*ones(N,1); % prior on the gamma3
gamma = kron(kron(gamma1,gamma2),gamma3);

keep_list = 1:N^3;
keep_list1 = 1:N;
keep_list2 = 1:N;% absolute position
keep_list3 = 1:N;% absolute position


norma = 1;
itra = 1;

r_absolute = 0;
r1 = 1e-2;
r2 = 1e-2;
r3 = 1e-2;
thres = 1e-4;

thres_inner = 1e-10;
R_max = 4000;
thres = 0;


while(itra < R_max + 1 && norma > thres) % do the iteration
    itra
    
    tic;

    gammal = zeros(N^3,1);
    gammal(keep_list) = gamma;
    gamma_old = gammal;
    

    % update the posterior mean and variance
    % compute the posterior
    [x_re,Sigma_x] = posterior_compute2(noise_var,gamma1,gamma2,gamma3,A1,A2,A3,A,y);
    
    Lambda = real(diag(Sigma_x)) + abs(x_re).^2;Sigma_x_diag = real(diag(Sigma_x));
    
    normi1 = 1;
    l1 = length(gamma1);
    l2 = length(gamma2);
    l3 = length(gamma3);

    while(normi1 > thres_inner)
        
        gamma_old_inner = kron(kron(gamma1,gamma2),gamma3);

        % update 1 and projection
        gamma1 = (l2*l3)^-1*(kron(eye(l1),kron(gamma2.^-1,gamma3.^-1)))'*Lambda;
        % update 2 and projection
        gamma2 = (l1*l3)^-1*(kron(gamma1.^-1,kron(eye(l2),gamma3.^-1)))'*Lambda;
        % update 3
        gamma3 = (l1*l2)^-1*(kron(gamma1.^-1,kron(gamma2.^-1,eye(l3))))'*Lambda;
        n1 = norm(gamma1);
        n2 = norm(gamma2);
        gamma1 = gamma1/n1; 
        gamma2 = gamma2/n2;
        gamma3 = gamma3*n1*n2;
        
        normi1 = norm(kron(kron(gamma1,gamma2),gamma3) - gamma_old_inner);        
    end
    
    time = time + toc;   

        
    gammal = zeros(N^3,1);
    gamma = kron(gamma1,kron(gamma2,gamma3));
    gammal(keep_list) = gamma;
    
    norma = norm(gammal - gamma_old)/norm(gamma_old);
    norm(gammal - gamma_old)/norm(gamma_old);
    
    
    % compute error for convergence tracking
    errora(itra) = (norm(x_re - x)/norm(x))^2;(norm(x_re - x)/norm(x))^2
    
    %
    itra = itra + 1;
end


metrics={'error',errora;
    'hyperparameter',[gamma1,gamma2,gamma3];
    };


end

