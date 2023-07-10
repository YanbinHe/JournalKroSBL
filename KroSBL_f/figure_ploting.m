%% some pre-defined values
% this file aims to reproduce the figures in the paper
run('Configuration.m');


% [all_themes, all_colors] = GetColors();
% all_colors = all_themes{1,1};
all_colors = [0, 70/255, 222/255;
              0.6350, 0.0780, 0.1840;
              0, 0.4470, 0.7410
              0.9290, 0.6940, 0.1250
              255/255,0,1;
              0.49 0.18 0.56;];


line_type_set{1} = '->';
line_type_set{2} = '-<';
line_type_set{3} = '-o';
line_type_set{4} = '-x';
line_type_set{5} = '-+';
line_type_set{6} = '-<';
line_type_set{7} = '-.>';
line_type_set{8} = '-.<';
line_type_set{9} = '-.o';
line_type_set{10} = '-.x';
line_type_set{11} = '-.+';
line_type_set{12} = '-.<';

legend_type_set{1} = '>';
legend_type_set{2} = '<';
legend_type_set{3} = 'o';
legend_type_set{4} = 'x';
legend_type_set{5} = '+';
legend_type_set{6} = '<';

algo_name{1} = 'cSBL';
algo_name{2} = 'AM-KroSBL';
algo_name{3} = 'SVD-KroSBL';
algo_name{4} = 'KOMP';
algo_name{5} = 'KSBL';
algo_name{6} = 'OLS';
num_alg = 5;
ratio = M1*M2*M3/N^3;

fontsizeman = 20;
opengl hardware
%%
% data process and figures plot
% take results.mat as input
% preprocessing
trials = 100;

% first index: algorithm; second index: 3 values (nmse/srr/time); third
% index: change with different conditions
resultSaggre = zeros(6,3,length(SNR));
resultMaggre = zeros(6,3,length(M1));
resultKaggre = zeros(6,3,length(K));



for t = 1:trials
    filename = ['./results_paper/noise_compare_', num2str(t),'.mat'];
    load(filename)
    for algo = 1:num_alg % for each algorithm
        metric = 1;
        for s = 1:length(SNR)
            resultSaggre(algo,metric,s) = resultSaggre(algo,metric,s) + resultS{s}{algo,1}{metric,2};
        end

        for m = 1:length(M1)
            resultMaggre(algo,metric,m) = resultMaggre(algo,metric,m) + resultM{m}{algo,1}{metric,2};
        end

        for k = 1:length(K)
            resultKaggre(algo,metric,k) = resultKaggre(algo,metric,k) + resultK{k}{algo,1}{metric,2};
        end

        metric = 3;
        for s = 1:length(SNR)
            resultSaggre(algo,metric-1,s) = resultSaggre(algo,metric-1,s) + resultS{s}{algo,1}{metric,2};
        end

        for m = 1:length(M1)
            resultMaggre(algo,metric-1,m) = resultMaggre(algo,metric-1,m) + resultM{m}{algo,1}{metric,2};
        end

        for k = 1:length(K)
            resultKaggre(algo,metric-1,k) = resultKaggre(algo,metric-1,k) + resultK{k}{algo,1}{metric,2};
        end

        metric = 4;
        for s = 1:length(SNR)
            resultSaggre(algo,metric-1,s) = resultSaggre(algo,metric-1,s) + resultS{s}{algo,1}{metric,2};
        end

        for m = 1:length(M1)
            resultMaggre(algo,metric-1,m) = resultMaggre(algo,metric-1,m) + resultM{m}{algo,1}{metric,2};
        end

        for k = 1:length(K)
            resultKaggre(algo,metric-1,k) = resultKaggre(algo,metric-1,k) + resultK{k}{algo,1}{metric,2};
        end
        
    end
end

resultSaggre = resultSaggre/trials;
resultMaggre = resultMaggre/trials;
resultKaggre = resultKaggre/trials;
%% fig 1/2 (a) NMSE/SRR vs SNR
figure
metric = 1;
for algo_index = 1:num_alg
    plot(SNR,reshape(resultSaggre(algo_index,metric,:),[6,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
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
xlabel('SNR (dB)')
ylabel('NMSE')
% axis square

h1 = plot(SNR,reshape(resultSaggre(1,metric,:),[6,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(SNR,reshape(resultSaggre(2,metric,:),[6,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(SNR,reshape(resultSaggre(3,metric,:),[6,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(SNR,reshape(resultSaggre(4,metric,:),[6,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(SNR,reshape(resultSaggre(5,metric,:),[6,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','northeast','Interpreter','LaTex')

figure
metric = 2;
for algo_index = 1:num_alg
    plot(SNR,reshape(resultSaggre(algo_index,metric,:),[6,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
% set(gca, 'yscale', 'log');
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
xlabel('SNR (dB)')
ylabel('SRR')
% axis square

h1 = plot(SNR,reshape(resultSaggre(1,metric,:),[6,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(SNR,reshape(resultSaggre(2,metric,:),[6,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(SNR,reshape(resultSaggre(3,metric,:),[6,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(SNR,reshape(resultSaggre(4,metric,:),[6,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(SNR,reshape(resultSaggre(5,metric,:),[6,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','northeast','Interpreter','LaTex')

%% plot figure 1/2 (b) NMSE/SRR vs measurements
figure
metric = 1;
for algo_index = 1:num_alg
    plot(ratio,reshape(resultMaggre(algo_index,metric,:),[5,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
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
xlabel('Under-sampling ratio')
ylabel('NMSE')
% axis square

h1 = plot(ratio,reshape(resultMaggre(1,metric,:),[5,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(ratio,reshape(resultMaggre(2,metric,:),[5,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(ratio,reshape(resultMaggre(3,metric,:),[5,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(ratio,reshape(resultMaggre(4,metric,:),[5,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(ratio,reshape(resultMaggre(5,metric,:),[5,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','northeast','Interpreter','LaTex')

figure
metric = 2;
for algo_index = 1:num_alg
    plot(ratio,reshape(resultMaggre(algo_index,metric,:),[5,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
% set(gca, 'yscale', 'log');
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
xlabel('Under-sampling ratio')
ylabel('SRR')
% axis square

h1 = plot(ratio,reshape(resultMaggre(1,metric,:),[5,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(ratio,reshape(resultMaggre(2,metric,:),[5,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(ratio,reshape(resultMaggre(3,metric,:),[5,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(ratio,reshape(resultMaggre(4,metric,:),[5,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(ratio,reshape(resultMaggre(5,metric,:),[5,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','northeast','Interpreter','LaTex')
%% figure 1/2 (c) NMSE/SRR vs sparsity level
figure
metric = 1;
for algo_index = 1:num_alg
    plot(K,reshape(resultKaggre(algo_index,metric,:),[5,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
ylim([1e-3,5]);
% axis square
h1 = plot(K,reshape(resultKaggre(1,metric,:),[5,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(K,reshape(resultKaggre(2,metric,:),[5,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(K,reshape(resultKaggre(3,metric,:),[5,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(K,reshape(resultKaggre(4,metric,:),[5,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(K,reshape(resultKaggre(5,metric,:),[5,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','northeast','Interpreter','LaTex')
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
xlabel('Sparsity level')
ylabel('NMSE')

% create a new pair of axes inside current figure
axes('position',[.55 .2 .35 .4])
box on % put box around new pair of axes
indexOfInterest = 3:5; % range of t near perturbation
for algo_index = [1,2,3,5]
    plot(K(indexOfInterest),reshape(resultKaggre(algo_index,metric,indexOfInterest),[3,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :));
    hold on
end
axis tight
grid on
set(gca, 'yscale', 'log');
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(gca,'YTick',[])
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
%%
figure
metric = 2;
for algo_index = 1:num_alg
    plot(K,reshape(resultKaggre(algo_index,metric,:),[5,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
% set(gca, 'yscale', 'log');
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
xlabel('Sparsity level')
ylabel('SRR')
% axis square

h1 = plot(K,reshape(resultKaggre(1,metric,:),[5,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(K,reshape(resultKaggre(2,metric,:),[5,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(K,reshape(resultKaggre(3,metric,:),[5,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(K,reshape(resultKaggre(4,metric,:),[5,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(K,reshape(resultKaggre(5,metric,:),[5,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','northeast','Interpreter','LaTex')

%% plot computation time

% k = 2, m = 2
% that is K = 3, M = 8
% using the results in resultS

timeAggre = zeros(6,5);

for t = 1:trials % for different independent runs
    filename = ['./results_paper/noise_compare_', num2str(t),'.mat'];
    load(filename)
    for algo = 1:num_alg % for each algorithm
        for k = 1:length(K)
            timeAggre(algo,k) = timeAggre(algo,k) + resultK{k}{algo,1}{4,2};
        end
    end
end

timeAggre = timeAggre/trials;

figure
for algo_index = 1:num_alg
    plot(K,timeAggre(algo_index,:),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end
fontsizeman = 15;
legend('boxoff')
title('(a)')
grid on
set(gca, 'yscale', 'log');
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
xlabel('SNR (dB)')
ylabel('Seconds (s)')
%%
% subplot(5,1,1)
% x = kron(x1,kron(x2,x3));
% plot(x);
% hold on
% subplot(5,1,2)
% plot(real(metrics_am{2,2}))
% hold on
% subplot(5,1,3)
% plot(metrics_csbl{2,2})
% hold on
% subplot(5,1,4)
% plot(metrics_omp{2,2})
% hold on
% subplot(5,1,5)
% plot(real(metrics_svd{2,2}))
%% Convergence performance
load('simu_results_con.mat');
figure
plot(metrics_am_con{1,2},'-','Color',all_colors(2, :),'Display',algo_name{2});
hold on
plot(metrics_svd_con{1,2},'-','Color',all_colors(3, :),'Display',algo_name{3});
hold on
plot(metrics_sota_con{1,2},'-','Color',all_colors(5, :),'Display',algo_name{5});
hold on
plot(metrics_csbl_con{1,2},'-','Color',all_colors(1, :),'Display',algo_name{1});
grid on

slot1 = 2.^[1:10];
slot2 = 2.^[1:8];
con_am = plot(slot1,metrics_am_con{1,2}(slot1),'<','Color',all_colors(2, :));
hold on
con_svd = plot(slot2,metrics_svd_con{1,2}(slot2),'o','Color',all_colors(3, :));
hold on
con_sota = plot(slot2,metrics_sota_con{1,2}(slot2),'+','Color',all_colors(5, :));
hold on
con_sbl = plot(slot2,metrics_csbl_con{1,2}(slot2),'>','Color',all_colors(1, :));
grid on

xlim([0,2500])
ylim([0.004,1])
yc = yline(metrics_csbl_con{1,2}(end),'--','NMSE = 0.0234','LineWidth',3);
hold on
ys = yline(metrics_svd_con{1,2}(end),'--','NMSE = 0.0066','LineWidth',3);
hold on
yso = yline(metrics_sota_con{1,2}(end),'--','NMSE = 0.0353','LineWidth',3);
hold on
ya = yline(metrics_am_con{1,2}(end),'--','NMSE = 0.0069','LineWidth',3);
ys.LabelVerticalAlignment = 'bottom';
yc.LabelHorizontalAlignment = 'left';
ys.LabelHorizontalAlignment = 'left';
ya.LabelHorizontalAlignment = 'left';
yso.LabelHorizontalAlignment = 'left';
yc.FontWeight = 'bold';
ys.FontWeight = 'bold';
ya.FontWeight = 'bold';
yso.FontWeight = 'bold';
yc.FontSize = 16;
ys.FontSize = 16;
ya.FontSize = 16;
yso.FontSize = 16;

legend([con_am,con_svd,con_sbl,con_sota],algo_name{2},algo_name{3},algo_name{1},algo_name{5},'Interpreter','LaTex');
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
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
xlabel('EM iteration number')
ylabel('NMSE')
legend('boxoff')
%% Convergence performance with prune
figure
plot(metrics_am_con2{1,2},'-','Color',all_colors(2, :),'Display',algo_name{2});
hold on
plot(metrics_svd_con2{1,2},'-','Color',all_colors(3, :),'Display',algo_name{3});
hold on
plot(metrics_sota_con2{1,2},'-','Color',all_colors(5, :),'Display',algo_name{5});
hold on
plot(metrics_csbl_con2{1,2},'-','Color',all_colors(1, :),'Display',algo_name{1});
grid on

slot1 = 2.^[1:10];
slot2 = 2.^[1:8];
con_am = plot(slot1,metrics_am_con2{1,2}(slot1),'<','Color',all_colors(2, :));
hold on
con_svd = plot(slot2,metrics_svd_con2{1,2}(slot2),'o','Color',all_colors(3, :));
hold on
con_sota = plot(slot2,metrics_sota_con2{1,2}(slot2),'+','Color',all_colors(5, :));
hold on
con_sbl = plot(slot2,metrics_csbl_con2{1,2}(slot2),'>','Color',all_colors(1, :));
grid on

xlim([0,500])
ylim([0.004,1])
yc = yline(metrics_csbl_con2{1,2}(end),'--','NMSE = 0.0234','LineWidth',3);
hold on
ys = yline(metrics_svd_con2{1,2}(end),'--','NMSE = 0.0066','LineWidth',3);
hold on
em = xline(150,'-','EM ends','LineWidth',3);
% yso = yline(metrics_sota_con2{1,2}(end),'--','NMSE = 0.0353','LineWidth',3);
% hold on
% ya = yline(metrics_am_con2{1,2}(end),'--','NMSE = 0.0069','LineWidth',3);
ys.LabelVerticalAlignment = 'bottom';
yc.LabelHorizontalAlignment = 'left';
ys.LabelHorizontalAlignment = 'left';
% ya.LabelHorizontalAlignment = 'left';
% yso.LabelHorizontalAlignment = 'left';
em.LabelHorizontalAlignment = 'left';
yc.FontWeight = 'bold';
ys.FontWeight = 'bold';
% ya.FontWeight = 'bold';
% yso.FontWeight = 'bold';
em.FontWeight = 'bold';
yc.FontSize = 16;
ys.FontSize = 16;
% ya.FontSize = 16;
% yso.FontSize = 16;
em.FontSize = 16;

legend([con_am,con_svd,con_sbl,con_sota],algo_name{2},algo_name{3},algo_name{1},algo_name{5},'Interpreter','LaTex');
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
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
xlabel('EM iteration number')
ylabel('NMSE')
legend('boxoff')
%% local minima

figure
local_minimum_am = kron(metrics_am_con{2,2}(:,1),kron(metrics_am_con{2,2}(:,2),metrics_am_con{2,2}(:,3)));
local_minimum_svd = kron(metrics_svd_con{2,2}(:,1),kron(metrics_svd_con{2,2}(:,2),metrics_svd_con{2,2}(:,3)));
subplot(3,1,1)

con_am = plot(local_minimum_am > 1e-4,'o','Color',all_colors(2, :),'Display',algo_name{2});
hold on
plot(abs(x) > 0,'kx');
grid on
title(algo_name{2})
xlim([1,6000])
set(gca, 'yscale', 'log');
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
xlabel('Position index')
set(gca,'yticklabel',[])
% ylabel('NMSE')
subplot(3,1,2)

con_svd = plot(local_minimum_svd> 1e-4,'o','Color',all_colors(3, :),'Display',algo_name{3});
hold on
plot(abs(x) > 0,'kx');
grid on
title(algo_name{3})
xlim([1,6000])
set(gca, 'yscale', 'log');
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
xlabel('Position index')
set(gca,'yticklabel',[])
% ylabel('NMSE')
subplot(3,1,3)

con_sbl = plot(metrics_csbl_con{2,2}> 1e-4,'o','Color',all_colors(1, :),'Display',algo_name{1});
hold on
posi = plot(abs(x) > 0,'kx');
grid on
title(algo_name{1})
xlim([1,6000])
set(gca, 'yscale', 'log');
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
xlabel('Position index')
set(gca,'yticklabel',[])
legend(posi,'True support');
% ylabel('NMSE')
%% rank plotting
load('simu_results.mat')
figure
plot(max(rank_svd),'Color',all_colors(3, :),'Display',algo_name{3});
hold on
plot(max(rank_am),'Color',all_colors(2, :),'Display',algo_name{2});
hold on
plot(max(rank_csbl),'Color',all_colors(1, :),'Display',algo_name{1});
hold on

plot_rank_svd = max(rank_svd);
plot_rank_am = max(rank_am);
plot_rank_csbl = max(rank_csbl);
h3 = plot(1:10:150,plot_rank_svd(1:10:150),legend_type_set{3},'Color',all_colors(3, :));
hold on
h2 = plot(1:10:150,plot_rank_am(1:10:150),legend_type_set{2},'Color',all_colors(2, :));
hold on
h1 = plot(1:10:150,plot_rank_csbl(1:10:150),legend_type_set{1},'Color',all_colors(1, :));
grid on
ylim([1,20]);
set(gca, 'yscale', 'log');
% set(gca, 'xscale', 'log');
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
xlabel('EM iteration number')
ylabel('Rank')

legend([h1 h2 h3],{algo_name{1},algo_name{2},algo_name{3}},'Location','northeast','Interpreter','LaTex')
