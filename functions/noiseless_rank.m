% noiseless rank
M1_local = 14;
K_local = 4;% sparsity level
R_max_local = 150;
       
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

y = kron(y1,kron(y2,y3));

% measurements
A1 = A1_ori;
A2 = A2_ori;
A3 = A3_ori;

A = kron(A1,kron(A2,A3));
%% classic SBL
[rank_csbl] = classicSBL_noiseless(y,A,N,R_max_local,x);
%% AM_KroSBL
[rank_am] = am_kroSBL_noiseless(y1,y2,y3,A1,A2,A3,N,R_max_local,x);
%% SVD_KroSBL
[rank_svd] = svd_kroSBL_noiseless(y1,y2,y3,A1,A2,A3,N,R_max_local,x);
