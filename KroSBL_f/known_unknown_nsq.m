clc
clear

addpath('./functions') 
%%
N = 15; % the dimension of sparse vector
M2 = 12; % the number of measurements per sparse vector
M3 = 15;
% noise convergence
M1_local = 14;
SNR_local = 10^(30/10); % in dB
K_local = 4;% sparsity level
R_max_local = 500;
func_ctrl = 1;
%% generate the sparse vector for testing
% generate measuring dictionaries
A1_ori = randn(M1_local,N);
A2_ori = randn(M2,N);
A3_ori = randn(M3,N);

supp1 = randsample(1:N, K_local);
supp2 = randsample(1:N, K_local);
supp3 = randsample(1:N, K_local);
%%
b1 = zeros(N,1);
b2 = zeros(N,1);
b3 = zeros(N,1);
b1(supp1) = ones(K_local,1);
b2(supp2) = ones(K_local,1);
b3(supp3) = ones(K_local,1);

b = kron(b1,kron(b2,b3));
suppTrue = find(abs(b)>0);

x = zeros(N^3,1);
x(suppTrue) = sqrt(0.05)*randn(K_local^3,1);

% measurements
A1 = A1_ori;
A2 = A2_ori;
A3 = A3_ori;

A = kron(A1,kron(A2,A3));

y_ori = A*x;

signal_power = norm(y_ori)^2/length(y_ori); %average signal power per asymbol

% SNR noise
noise_var = (signal_power)/SNR_local;
noise = sqrt(noise_var)*randn(size(y_ori));
y = y_ori + noise;
%% unknown noise variance
%% AM_KroSBL
[metrics_am_con] = am_kroSBL_con_nsq(y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% SVD_KroSBL
[metrics_svd_con] = svd_kroSBL_con_nsq(y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% KSBL
[metrics_sota_con] = sota_kroSBL_con_nsq(y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% classic SBL
[metrics_csbl_con] = classicSBL_con_nsq(y,A,N,R_max_local,x,func_ctrl);
%% known noise variance
%% AM_KroSBL
[metrics_am_con3] = am_kroSBL_con3(noise_var,y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% SVD_KroSBL
[metrics_svd_con3] = svd_kroSBL_con3(noise_var,y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% KSBL
[metrics_sota_con3] = sota_kroSBL_con3(noise_var,y,A1,A2,A3,A,N,R_max_local,x,func_ctrl);
%% classic SBL
[metrics_csbl_con3] = classicSBL_con3(noise_var,y,A,N,R_max_local,x,func_ctrl);