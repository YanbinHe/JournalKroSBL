function [metrics] = classicSBL_con2(noise_var,y,A,N,R_max,x,func_ctrl)

time = 0;
A_ori = A;
gammac = 1*ones(N^3,1); % initialize the prior
x_re = 1*ones(N^3,1);

keep_listc = 1:N^3;

normc = 1;
itrc = 1;

r_rele = 5*1e-4;
thres = 1e-4;

if func_ctrl.noisy_convergence
    R_max = 500;
    thres = 0;
end

while(itrc < R_max + 1 && normc > thres) % do the iteration
    itrc
    
    tic;
    
    gammal = zeros(N^3,1);
    gammal(keep_listc) = gammac;
    gammac_old = gammal;
    
        if func_ctrl.noisy_convergence

            if min(abs(gammac)) < r_rele*max(abs(gammac))% 
                indexc = find(gammac > r_rele*max(abs(gammac))); % should keep 
                gammac = gammac(indexc);
                A = A(:,indexc);
                keep_listc = keep_listc(indexc);   
            end
    
        end

    % update the posterior
    Gammac = diag(gammac);
    invPhic = A'* (noise_var * eye(size(A,1)) + A * Gammac * A')^-1 * A;
    Sigma_x = Gammac - Gammac * invPhic * Gammac;
    x_re = noise_var^-1 * Sigma_x * A' * y;
    
    % updating gamma
    Sigma_x_diag = real(diag(Sigma_x));
    gammac = Sigma_x_diag + abs(x_re).^2;

    time = time + toc;
    
    gammal = zeros(N^3,1);
    gammal(keep_listc) = gammac;

    normc = norm(gammal-gammac_old)/norm(gammac_old);
    norm(gammal-gammac_old)/norm(gammac_old)
    
    % compute error for convergence tracking
    
    x_rel = zeros(N^3,1);
    x_rel(keep_listc) = x_re;

    errorc(itrc) = norm(x_rel - x)/norm(x);
    
    itrc = itrc + 1; 
end



if ~func_ctrl.noisy_convergence
    srr = recover_rate(x_rel,x);
    metrics={'error',errorc;
             'vector',x_rel;
             'support recovery rate',srr;
             'time',time
             };
else
    metrics={'error',errorc;
        'hyperparameter',gammal;
        };
end

end

