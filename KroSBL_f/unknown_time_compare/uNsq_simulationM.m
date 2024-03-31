% This file is fed with different s/m/k values to conduct simulation
% not a function

% measurements
A1 = A1_ori(:,:);
A2 = A2_ori(:,:);
A3 = A3_ori(:,:);

y_noiseless = y_ori;
A = kron(A1,kron(A2,A3));


% SNR noise
noise_var = (signal_power)/SNR_10(s);
noise = sqrt(noise_var)*randn(size(y_noiseless));
y = y_noiseless + noise;

% implementation of each algorithm
%% classic SBL
[metrics_csbl] = classicSBL_un(y,A,N,R_max,x,func_ctrl);
%% AM_KroSBL
[metrics_am] = am_kroSBL_un(y,A1,A2,A3,A,N,R_max,x,func_ctrl);
%% SVD_KroSBL
[metrics_svd] = svd_kroSBL_un(y,A1,A2,A3,A,N,R_max,x,func_ctrl);
%% KroSBL_sota
[metrics_sota] = sota_kroSBL_un(y,A1,A2,A3,A,N,R_max,x,func_ctrl);
%% data saving
result = cell(6,1);
result{1} = metrics_csbl;
result{2} = metrics_am;
result{3} = metrics_svd;
result{5} = metrics_sota;
resultM{m} = result;