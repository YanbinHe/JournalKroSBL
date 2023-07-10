% symbol error rate check
clc
clear

%% simluation on the SBL based IRS channel estimation with spread
% predefined values
ms_ante = 4; % the number of antennas at user equippment side
bs_ante = 16; % the number of antennas at bs
irs_ele = 16^2; % the number of irs elements
Res1 = 18; % the number of atoms in the dicitionay
% trans_power = 30; % in dBm
% trans_power_10 = 10^(trans_power/10 - 3); % in Watt

overheadp = [6,10,14,18]; % the number of overhead w.r.t pilot signals
overheadi = [6,10,14,18]; % the number of overhead w.r.t reflecting pattern

angle_spread = 3; % the spread of angles in the domain of cos
angle = linspace(-1,1-2/Res1,Res1); 
% d_msirs = 10; % distance between user and irs
% d_irsbs = 30; % distance between irs and base station

SNRl = [-5,0,5,10,15,20,25,30]; % in dB
SNRl_10 = 10.^(SNRl/10);

thres = 1e-4;
global thres_inner;
thres_inner = 1e-3;

numItr = 150; % the number of alternative optimization iterations
klevel = 30;
%% simulation settings
% X = 1/sqrt(2*ms_ante)*(randn(ms_ante,max(overheadp))+1i*randn(ms_ante,max(overheadp))); % pilot signal
irs_pattern = 1/sqrt(irs_ele)*((rand(irs_ele,max(overheadi)) > 0.5)*2-1); 
%%
% different channel realization, stay constant within the coherence time
Smpl = randsample([0:Res1-angle_spread], 1);
AoD_ms = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoA_irs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoD_irs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoA_bs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;  

% generate the true support, for the support recovery computation
g2 = zeros(Res1,1);
g2(AoA_bs) = 1;

g1 = zeros(Res1,1);
g1(AoD_ms) = 1;

gi = zeros(Res1,1);
indexi = [(AoA_irs(1)-AoD_irs(end)):1:(AoA_irs(end)-AoD_irs(1))] + Res1;
indexmodi = find(indexi > Res1);
indexi(indexmodi) = indexi(indexmodi) - Res1;
gi(indexi) = 1; 
gi = (fliplr(gi'))';

suppTrue = kron(gi,kron(g1,g2));
%%
% array response vector generate (matrix)
ar_AoD_ms = generate_steering(ms_ante,angle(AoD_ms));
ar_AoA_irs = generate_steering(irs_ele,angle(AoA_irs));
ar_AoD_irs = generate_steering(irs_ele,angle(AoD_irs));
ar_AoA_bs = generate_steering(bs_ante,angle(AoA_bs));

% path gain CN(0,1)
alpha_1 = (1*randn(1,angle_spread) + 1i*randn(1,angle_spread))/sqrt(2); % from ms to irs
alpha_2 = (1*randn(1,angle_spread) + 1i*randn(1,angle_spread))/sqrt(2); % from irs to bs

H2 = sqrt(bs_ante*irs_ele/angle_spread)*ar_AoA_bs *diag(alpha_2)*ar_AoD_irs';
H1 = sqrt(ms_ante*irs_ele/angle_spread)*ar_AoA_irs*diag(alpha_1)*ar_AoD_ms';

% make sure that somehow the irs should point to the receiver, otherwise
% the received signal is too weak. But this is only for the simulation.
% This is not the case when it comes to real world.
for i = 1:max(overheadi)
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    while norm(Htrue(:,:,i),'fro')<1e-5
        irs_pattern(:,i) = 1/sqrt(irs_ele)*(rand(irs_ele,1) > 0.5)*2-1; 
        Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    end
end
%%
% SNR and overhead settings
K1 = overheadp(1);
K2 = overheadi(1);

%%
re = [];
ser = 0;
SNR = 35;
nv = 10^(-SNR/10);

for i = 1:K2
    Htrue = H2*diag(irs_pattern(:,1))*H1;
    noise = sqrt(nv / 2) * (randn(bs_ante,ms_ante) + 1i*randn(bs_ante,ms_ante));
    H = Htrue + noise;
    norm(H - Htrue,'fro')/norm(Htrue,'fro')
    
    [s,sp] = ser_compute(H,Htrue,SNR,1e4);
    
    ser(i) = s;
    spa(i) = sp;
end

re = sum(ser)/K2


% plot([-5:5:100],ser)



