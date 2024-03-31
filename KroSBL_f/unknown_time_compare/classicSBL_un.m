function [metrics] = classicSBL_un(y,A,N,R_max,x,func_ctrl)

time = 0;
A_ori = A;
gammac = 1*ones(N^3,1); % initialize the prior
num_M = length(y);

keep_listc = 1:N^3;
noise_var = 1e-1*var(y);

normc = 1;
itrc = 1;

r_rele = 1e-3;
thres = 1e-4;

while(itrc < R_max + 1 && normc > thres) % do the iteration
    itrc;
    
    tic;
    
    gammal = zeros(N^3,1);
    gammal(keep_listc) = gammac;
    gammac_old = gammal;
    
        if ~func_ctrl.noisy_convergence

            if min(abs(gammac)) < r_rele*max(abs(gammac))% 
                indexc = find(gammac > r_rele*max(abs(gammac))); % should keep 
                gammac = gammac(indexc);
                A = A(:,indexc);
                keep_listc = keep_listc(indexc);   
            end
    
        end

        gammad = gammac;
    % update the posterior
    Gammac = diag(gammac);
    invPhic = A'* pinv(noise_var * eye(size(A,1)) + A * Gammac * A') * A;
    Sigma_x = Gammac - Gammac * invPhic * Gammac;
    x_re = noise_var^-1 * Sigma_x * A' * y;
    
    % updating gamma
    Sigma_x_diag = real(diag(Sigma_x));
    gammac = Sigma_x_diag + abs(x_re).^2;

    % update noise variance
    temp = sum(1 - Sigma_x_diag./gammad);
    noise_var = (norm(y - A*x_re)^2 + noise_var*temp)/num_M;

    time = time + toc;
    
    gammal = zeros(N^3,1);
    gammal(keep_listc) = gammac;

    normc = norm(gammal-gammac_old)/norm(gammac_old);
    norm(gammal-gammac_old)/norm(gammac_old);
    itrc = itrc + 1; 
end

% re-estimate
Gammac = diag(gammal);
invPhic = A_ori'* (noise_var * eye(size(A_ori,1)) + A_ori * Gammac * A_ori')^-1 * A_ori;
Sigma_x = Gammac - Gammac * invPhic * Gammac;
x_rel = noise_var^-1 * Sigma_x * A_ori' * y;

errorc = (norm(x_rel - x)/norm(x))^2

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

