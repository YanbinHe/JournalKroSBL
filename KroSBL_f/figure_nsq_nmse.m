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

%% plot
ratio = M1*M2*M3/N^3;

trials = 200;
lenM = length(M1);
errorAggre = zeros(6,lenM);
srrAggre = zeros(6,lenM);
timeAggre = zeros(6,lenM);
for t = 1:trials % for different independent runs
    filename = ['./resultsa/unknown_nsq_nmse_m_', num2str(t),'.mat'];
    load(filename)
    for m = 1:lenM
        for algo = [1,2,3,5] % for each algorithm
            errorAggre(algo,m) = errorAggre(algo,m) + resultM{m}{algo,1}{1,2};
            srrAggre(algo,m) = srrAggre(algo,m) + resultM{m}{algo,1}{3,2};
            timeAggre(algo,m) = timeAggre(algo,m) + resultM{m}{algo,1}{4,2};
        end
    end
end

errorAggre = errorAggre/trials;
srrAggre = srrAggre/trials;
timeAggre = timeAggre/trials;

figure
for algo_index = [1,2,3,5]
    plot(ratio,errorAggre(algo_index,:),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

legend('boxoff')
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
xlabel('Under-sampling ratio','Interpreter','latex')
ylabel('NMSE','Interpreter','latex')

h1 = plot(ratio,errorAggre(1,:),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(ratio,errorAggre(2,:),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(ratio,errorAggre(3,:),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h5 = plot(ratio,errorAggre(5,:),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{5}},'Location','northeast','Interpreter','LaTex')

figure
for algo_index = [1,2,3,5]
    plot(ratio,srrAggre(algo_index,:),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

legend('boxoff')
grid on
% set(gca, 'yscale', 'log');
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
xlabel('Under-sampling ratio','Interpreter','latex')
ylabel('SRR','Interpreter','latex')
ylim([0,1])

figure
for algo_index = 1:num_alg
    plot(ratio,timeAggre(algo_index,:),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
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
xlabel('SNR (dB)','Interpreter','LaTex')
ylabel('Seconds (s)','Interpreter','LaTex')