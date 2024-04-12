% Default values for configurations
% we fix the dimension of the second/third matrix
addpath('./functions') 
%%


N = 15; % the dimension of sparse vector
M1 = [6,8,10,12,14];
M2 = 12; % the number of measurements per sparse vector
M3 = 15;
K = [2,3,4,5,6]; % the number of non-zero entries
SNR = [5,10,15,20,25,30]; % SNR values
SNR_10 = 10.^(SNR/10);
AVG = 200; % the number of trials

R_max = 150; % the maximum number of iterations for SBL/KroSBL
result = [];
