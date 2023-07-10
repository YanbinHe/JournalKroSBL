% noise convergence
M1_local = 14;
SNR_local = 10^(30/10); % in dB
K_local = 4;% sparsity level
R_max_local = 500;
       
%% generate the sparse vector for testing
% generate measuring dictionaries
A1_ori = randn(M1_local,N);
A2_ori = randn(M2,N);
A3_ori = randn(M3,N);

x1 = zeros(N,1);
x2 = zeros(N,1);
x3 = zeros(N,1);

supp1 = randsample(1:N, K_local);
x1(supp1) = randn(K_local,1);
supp2 = randsample(1:N, K_local);
x2(supp2) = randn(K_local,1);
supp3 = randsample(1:N, K_local);
x3(supp3) = randn(K_local,1);

x = kron(x1,kron(x2,x3));
suppTrue = find(abs(x)>0);

y1 = A1_ori*x1;
y2 = A2_ori*x2;
y3 = A3_ori*x3;

y_ori = kron(y1,kron(y2,y3));

signal_power = norm(y_ori)^2/length(y_ori); %average signal power per asymbol

% measurements
A1 = A1_ori;
A2 = A2_ori;
A3 = A3_ori;

A = kron(A1,kron(A2,A3));

% SNR noise
noise_var = (signal_power)/SNR_local;
noise = sqrt(noise_var)*randn(size(y_ori));
y = y_ori + noise;
%% AM_KroSBL
[metrics_am_con] = am_kroSBL_con(noise_var,y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% SVD_KroSBL
[metrics_svd_con] = svd_kroSBL_con(noise_var,y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% KSBL
[metrics_sota_con] = sota_kroSBL_con(noise_var,y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% classic SBL
[metrics_csbl_con] = classicSBL_con(noise_var,y,A,N,R_max_local,x,func_ctrl);