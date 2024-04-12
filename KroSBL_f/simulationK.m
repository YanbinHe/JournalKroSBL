% This file is fed with different s/m/k values to conduct simulation
% not a function

k = 1;
supp1 = randsample(1:N, K(k));
supp2 = randsample(1:N, K(k));
supp3 = randsample(1:N, K(k));

x1 = zeros(N,1);
x2 = zeros(N,1);
x3 = zeros(N,1);

for k = 1:lenK

    if k == 1

        x1(supp1) = randn(K(k),1);
        x2(supp2) = randn(K(k),1);
        x3(supp3) = randn(K(k),1);

    else

        zo1 = find(x1 == 0);
        zo2 = find(x2 == 0);
        zo3 = find(x3 == 0);

        new1 = randsample(zo1,1);
        new2 = randsample(zo2,1);
        new3 = randsample(zo3,1);

        x1(new1) = randn();
        x2(new2) = randn();
        x3(new3) = randn();
    end

    x = kron(x1,kron(x2,x3));
    suppTrue = find(abs(x)>0);

    y1 = A1_ori(1:M1(m),:)*x1;
    y2 = A2_ori*x2;
    y3 = A3_ori*x3;

    y_ori = kron(y1,kron(y2,y3));

    signal_power = norm(y_ori)^2/length(y_ori); %average signal power per asymbol

    % measurements
    A1 = A1_ori(1:M1(m),:);
    A2 = A2_ori(:,:);
    A3 = A3_ori(:,:);

    A = kron(A1,kron(A2,A3));

    % SNR noise
    noise_var = (signal_power)/SNR_10(s);
    noise = sqrt(noise_var)*randn(size(y_ori));
    y = y_ori + noise;

    % implementation of each algorithm
    %% classic SBL
    [metrics_csbl] = classicSBL(noise_var,y,A,N,R_max,x,func_ctrl);
    %% AM_KroSBL
    [metrics_am] = am_kroSBL(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl);
    %% SVD_KroSBL
    [metrics_svd] = svd_kroSBL(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl);
    %% KroSBL_sota
    [metrics_sota] = sota_kroSBL(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl);
    %% Kronecker-OMP
    [metrics_omp] = KroOMP(A1,A2,A3,y,OMP_level,epsilon,x);
    %% benchmark-oracle LS
    [metrics_ols] = OLS(y,A,N,suppTrue,x);
    %% data saving
    result = cell(6,1);
    result{1} = metrics_csbl;
    result{2} = metrics_am;
    result{3} = metrics_svd;
    result{4} = metrics_omp;
    result{5} = metrics_sota;
    result{6} = metrics_ols;
    resultK{k} = result;
end
% % This file is fed with different s/m/k values to conduct simulation
% % not a function
% 
% x1 = zeros(N,1);
% x2 = zeros(N,1);
% x3 = zeros(N,1);
% 
% supp1 = randsample(1:N, K(k));
% x1(supp1) = randn(K(k),1);
% supp2 = randsample(1:N, K(k));
% x2(supp2) = randn(K(k),1);
% supp3 = randsample(1:N, K(k));
% x3(supp3) = randn(K(k),1);
% 
% x = kron(x1,kron(x2,x3));
% suppTrue = find(abs(x)>0);
% 
% y1 = A1_ori(1:M1(m),:)*x1;
% y2 = A2_ori*x2;
% y3 = A3_ori*x3;
% 
% y_ori = kron(y1,kron(y2,y3));
% 
% signal_power = norm(y_ori)^2/length(y_ori); %average signal power per asymbol
% 
% % measurements
% A1 = A1_ori(1:M1(m),:);
% A2 = A2_ori(:,:);
% A3 = A3_ori(:,:);
% 
% A = kron(A1,kron(A2,A3));
% 
% % SNR noise
% noise_var = (signal_power)/SNR_10(s);
% noise = sqrt(noise_var)*randn(size(y_ori));
% y = y_ori + noise;
% 
% % implementation of each algorithm
% %% classic SBL
% [metrics_csbl] = classicSBL(noise_var,y,A,N,R_max,x,func_ctrl);
% %% AM_KroSBL
% [metrics_am] = am_kroSBL(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl);
% %% SVD_KroSBL
% [metrics_svd] = svd_kroSBL(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl);
% %% KroSBL_sota
% [metrics_sota] = sota_kroSBL(noise_var,y,A1,A2,A3,A,N,R_max,x,func_ctrl);
% %% Kronecker-OMP
% [metrics_omp] = KroOMP(A1,A2,A3,y,OMP_level,epsilon,x);
% %% benchmark-oracle LS
% [metrics_ols] = OLS(y,A,N,suppTrue,x);
% %% data saving
% result = cell(6,1);
% result{1} = metrics_csbl;
% result{2} = metrics_am;
% result{3} = metrics_svd;
% result{4} = metrics_omp;
% result{5} = metrics_sota;
% result{6} = metrics_ols;
% resultK{k} = result;