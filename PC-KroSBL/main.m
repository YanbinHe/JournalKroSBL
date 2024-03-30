clc
clear
addpath('functions')
%% simluation on the SBL based IRS channel estimation with angle spread
% predefined values
ms_ante = 6; % the number of antennas at user equippment side
bs_ante = 16; % the number of antennas at bs
irs_ele = 16^2; % the number of irs elements
Res1 = 18; % the number of atoms in the dicitionay
timeslots = 1e5;
phase_offset = 0;

overheadp = 3; % the number of overhead w.r.t pilot signals
overheadi = 6; % the number of overhead w.r.t reflecting pattern

angle_spread = 3; % the spread of angles in the domain of cos
angle = linspace(-1,1-2/Res1,Res1);

SNRl = [4,8,12,16,20]; % in dB
SNRl_10 = 10.^(SNRl/10);

thres = 1e-4;
global thres_inner;
thres_inner = 1e-4;

numItr = 150; % the number of alternative optimization iterations
klevel = 50;

symbolmatrix = pskmod([0:4-1], 4, phase_offset);
%% simulation settings

JJ = 5;
AVG = 50;% the number of trials

X = 1/sqrt(ms_ante)*pilot_gen(ms_ante,max(overheadp));
irs_pattern = 1/sqrt(irs_ele)*((rand(irs_ele,max(overheadi)) > 0.5)*2-1);

%%
for avg = 1 : AVG % for AVG trials
    % different channel realization, stay constant within the coherence time
    error = [];
    time = [];
    ser = [];

    AoD_ms = randsample([1:Res1], 1);
    g1 = zeros(Res1,1);
    g1(AoD_ms) = 1;
    
    Smpl = randsample([0:Res1-angle_spread], 1);
    AoA_irs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;
    ga = zeros(Res1,1);
    ga(AoA_irs) = 1;
    
    AoD_irs = randsample([1:Res1], 1);
    gd = zeros(Res1,1);
    gd(AoD_irs) = 1;
    
    Smpl = randsample([0:Res1-angle_spread], 1);
    AoA_bs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;
    g2 = zeros(Res1,1);
    g2(AoA_bs) = 1;
    
    % generate the true support, for the support recovery computation
    gi = zeros(Res1,1);
    indexi = [(AoA_irs(1)-AoD_irs(end)):1:(AoA_irs(end)-AoD_irs(1))] + Res1;
    indexmodi = find(indexi > Res1);
    indexi(indexmodi) = indexi(indexmodi) - Res1;
    gi(indexi) = 1;
    gi = (fliplr(gi'))';
    
    suppTrue = kron(gi,kron(g1,g2));
    %%
    % dictionary generation
    A2 = generate_dict(bs_ante,Res1);
    A_irs_d = generate_dict(irs_ele,Res1);
    A_irs_a = generate_dict(irs_ele,Res1);
    A1 = generate_dict(ms_ante,Res1);
    
    % path gain CN(0,1)
    alpha_aco = (1*randn(angle_spread,1) + 1i*randn(angle_spread,1))/sqrt(2); % from ms to irs
    alpha_a = zeros(Res1,1);
    alpha_a(AoA_irs) = alpha_aco;
    
    alpha_2co = (1*randn(angle_spread,1) + 1i*randn(angle_spread,1))/sqrt(2); % from irs to bs
    alpha_2 = zeros(Res1,1);
    alpha_2(AoA_bs) = alpha_2co;
    
    H2 = sqrt(bs_ante*irs_ele/angle_spread)*A2*(g2.*alpha_2)*gd'*A_irs_d';
    H1 = sqrt(ms_ante*irs_ele/angle_spread)*A_irs_a*(ga.*alpha_a)*g1'*A1';
    
    % make sure that somehow the irs should point to the receiver, otherwise
    % the received signal is too weak. But this is only for the simulation.
    % This is not the case when it comes to real world.
    for i = 1:max(overheadi)
        Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
        while norm(Htrue(:,:,i),'fro')<1e-5
            irs_pattern(:,i) =  1/sqrt(irs_ele)*(rand(irs_ele,1) > 0.5)*2-1;
            Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
        end
    end
    
    % generate signal
    y_bar = kr(X.'*H1.',H2)*irs_pattern;
    signal_power = norm(vec(y_bar))^2/length(vec(y_bar)); %average signal power per asymbol
    %% start simulation
    for jj = 1:JJ % snr
        
        SNR_10 = SNRl_10(jj);
        noise_var = (signal_power)/SNR_10;
        
        
        % SNR and overhead settings
        K1 = overheadp;
        K2 = overheadi;
        
        x_pilot = X(:,1:K1);
        IRS = irs_pattern(:,1:K2);
        % different irs reflection patterns, the number is equal to #overheads,
        H_p1 = IRS.'*kr(A_irs_a.',A_irs_d').';
        H_p1 = H_p1(:,1:Res1);
        H_p1_ori = H_p1;
        H_p2 = x_pilot.'*conj(A1);
        H_p2_ori = H_p2;
        H = kron(kron(H_p1,H_p2),A2);
        %% SNR part
        % received noisy signals
        y_tilde = vec(y_bar(1:K1*bs_ante,1:K2));
        % generate noise
        noise = sqrt(noise_var / 2)*(randn(size(y_tilde))+1i*randn(size(y_tilde)));
        y = y_tilde + noise;
        %% part 2: channel estimation with different techniques and compute the symbol error rate
        % true channel for each IRS pattern
        for i = 1:K2
            Htrue(:,:,i) = vec(H2*diag(IRS(:,i))*H1);
        end
        
        inPut = randi([0, 4-1], 1, timeslots);
        Xun = pskmod(inPut, 4, phase_offset);
        
        %% different techniques
        %% classicSBL
        [error_csbl,time_csbl,~,g_csbl,Hre_csbl] = classicSBL(numItr,H,Res1,noise_var,y,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2,SNRl(jj),thres);
        % metrics compute
        error(7,jj) = error_csbl;
        time(7,jj) = time_csbl;
        ser(7,jj) = ser_compute(Hre_csbl,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
        %% PCSBL
        [error_pcsbl,time_pcsbl,~,g_pcsbl,Hre_pcsbl] = PCSBL(numItr,H,Res1,noise_var,y,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2,SNRl(jj),thres);
        % metrics compute
        error(4,jj) = error_pcsbl;
        time(4,jj) = time_pcsbl;
        ser(4,jj) = ser_compute(Hre_pcsbl,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
        %% SVD
        [error_svdsbl,time_svdsbl,~,gl_svdsbl,Hre_svdsbl] = kroSBL_svd(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2);
        % metrics compute
        error(1,jj) = error_svdsbl;
        time(1,jj) = time_svdsbl;
        ser(1,jj) = ser_compute(Hre_svdsbl,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
        %% Alternating
        [error_am,time_am,~,gl_am,Hre_am] = kroSBL_am(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2);
        error(2,jj) = error_am;
        time(2,jj) = time_am;
        ser(2,jj) = ser_compute(Hre_am,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
        %% PC alpha = 0.9
        [error_pc9,time_pc9,~,gl_pc9,Hre_pc9] = kroSBL_am_pc(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2,0.8,0.8);
        % metrics compute
        error(5,jj) = error_pc9;
        time(5,jj) = time_pc9;
        ser(5,jj) = ser_compute(Hre_pc9,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
        %% PC alpha = 0.5
        [error_pc5,time_pc5,~,gl_pc5,Hre_pc5] = kroSBL_am_pc(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2,0.5,0.5);
        % metrics compute
        error(6,jj) = error_pc5;
        time(6,jj) = time_pc5;
        ser(6,jj) = ser_compute(Hre_pc5,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
        
    end
    filename = ['./Results/compare_noisy_try_', num2str(avg),'.mat'];
    save(filename, 'error', 'time', 'ser')
end
%% Plot part
%% ii-overhead jj-snr preprocessing
trial = 50;
lengt = 5;
errorAggr = zeros(7,lengt); % 8 algorithms
timeAggr = zeros(7,lengt);
serAggr = zeros(7,lengt);
for t = 1:trial
    filename = ['./Results/compare_noisy_', num2str(t),'.mat'];
    load(filename)
    for algo_index = [1,2,4,5,6,7]
        errorAggr(algo_index,:) = errorAggr(algo_index,:) + error(algo_index,:);
        timeAggr(algo_index,:) = timeAggr(algo_index,:) + time(algo_index,:);
        serAggr(algo_index,:) = serAggr(algo_index,:) + ser(algo_index,:);
    end
end
errorAggr = errorAggr/trial;
timeAggr = timeAggr/trial;
serAggr = serAggr/trial;
%
fontsizeman = 20;
num_alg = 8;
          
all_colors = [0, 0.4470, 0.7410;
              0.6350, 0.0780, 0.1840;
              0 0 0;
              0.49 0.18 0.56
              0.9290, 0.6940, 0.1250
              255/255,0,1;
              0, 70/255, 222/255];


line_type_set{1} = '-o';
line_type_set{2} = '-<';
line_type_set{3} = '-+';
line_type_set{4} = '-*';
line_type_set{5} = '-^';
line_type_set{6} = '-v';
line_type_set{7} = '->';
line_type_set{8} = '-x';

legend_type_set{1} = 'o';
legend_type_set{2} = '<';
legend_type_set{3} = '+';
legend_type_set{4} = '*';
legend_type_set{5} = '^';
legend_type_set{6} = 'v';
legend_type_set{7} = '>';
legend_type_set{8} = 'x';


algo_name{1} = 'SVD-KroSBL';
algo_name{2} = 'AM-KroSBL';
algo_name{3} = 'KSBL';
algo_name{4} = 'PC-SBL';
algo_name{5} = 'PC-KroSBL $\beta$=0.8';
algo_name{6} = 'PC-KroSBL $\beta$=0.5';
algo_name{7} = 'cSBL';
algo_name{8} = 'OMP';
%
figure
for algo_index = [1,2,4,5,6,7]
    plot(SNRl,errorAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
    hold on
end

h1 = plot(SNRl,errorAggr(1,:),legend_type_set{1},'Display',algo_name{1},'Color',all_colors(1, :));
hold on
h2 = plot(SNRl,errorAggr(2,:),legend_type_set{2},'Display',algo_name{2},'Color',all_colors(2, :));
hold on
h4 = plot(SNRl,errorAggr(4,:),legend_type_set{4},'Display',algo_name{4},'Color',all_colors(4, :));
hold on
h5 = plot(SNRl,errorAggr(5,:),legend_type_set{5},'Display',algo_name{5},'Color',all_colors(5, :));
hold on
h6 = plot(SNRl,errorAggr(6,:),legend_type_set{6},'Display',algo_name{6},'Color',all_colors(6, :));
hold on
h7 = plot(SNRl,errorAggr(7,:),legend_type_set{7},'Display',algo_name{7},'Color',all_colors(7, :));
hold on
legend([h1 h2 h4 h5 h6 h7],{algo_name{1},algo_name{2},algo_name{4},algo_name{5},algo_name{6},algo_name{7}},'Location','northeast','Interpreter','LaTex')

grid on
set(gca, 'yscale', 'log');
% set(gca, 'xscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlim([4,20])
xticks(SNRl)
xlabel('SNR (dB)','Interpreter','LaTex')
ylabel('NMSE','Interpreter','LaTex')

%
figure
for algo_index = [1,2,4,5,6,7]
    plot(SNRl,serAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
    hold on
end

grid on
set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlim([4,20])
xticks(SNRl)
xlabel('SNR (dB)','Interpreter','LaTex')
ylabel('SER','Interpreter','LaTex')

figure
for algo_index = [1,2,4,5,6,7]
    plot(SNRl,timeAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
    hold on
end

grid on
set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlim([4,20])
xticks(SNRl)
xlabel('SNR (dB)')
ylabel('NMSE')