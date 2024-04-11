%% Convergence performance
load('simu_nsq_vs_nonsq.mat')
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
fontsizeman = 20;
figure

plot(metrics_am_con{1,2},'-','Color',all_colors(2, :),'Display',algo_name{2});
hold on
plot(metrics_svd_con{1,2},'-','Color',all_colors(3, :),'Display',algo_name{3});
hold on
plot(metrics_sota_con{1,2},'-','Color',all_colors(5, :),'Display',algo_name{5});
hold on
plot(metrics_csbl_con{1,2},'-','Color',all_colors(1, :),'Display',algo_name{1});
grid on
hold on

slot1 = 2.^[0:10];
slot2 = 2.^[0:8];
slot3 = 2.^[0:7];
con_am = plot(slot1,metrics_am_con{1,2}(slot1),'<','Color',all_colors(2, :));
hold on
con_svd = plot(slot3,metrics_svd_con{1,2}(slot3),'o','Color',all_colors(3, :));
hold on
con_sota = plot(slot2,metrics_sota_con{1,2}(slot2),'+','Color',all_colors(5, :));
hold on
con_sbl = plot(slot2,metrics_csbl_con{1,2}(slot2),'>','Color',all_colors(1, :));
grid on
hold on

xlim([0,4000])
ylim([0.00001,10])
yc = yline(metrics_csbl_con{1,2}(end),'--','NMSE = 0.0016','LineWidth',3);
hold on
ys = yline(metrics_svd_con{1,2}(end),'--','NMSE = 0.000039','LineWidth',3);
hold on
yso = yline(metrics_sota_con{1,2}(end),'--','NMSE = 0.00037','LineWidth',3);
hold on
ya = yline(metrics_am_con{1,2}(end),'--','NMSE = 0.000038','LineWidth',3);
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
title('Unknown noise variance')
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
% Convergence performance
% load('simu_results_con.mat');
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
fontsizeman = 20;
figure

plot(metrics_am_con3{1,2},':','Color',all_colors(2, :),'Display',algo_name{2});
hold on
plot(metrics_svd_con3{1,2},':','Color',all_colors(3, :),'Display',algo_name{3});
hold on
plot(metrics_sota_con3{1,2},':','Color',all_colors(5, :),'Display',algo_name{5});
hold on
plot(metrics_csbl_con3{1,2},':','Color',all_colors(1, :),'Display',algo_name{1});
grid on
hold on
slot1 = 2.^[0:10];
slot2 = 2.^[0:8];
slot3 = 2.^[0:7];
con_am = plot(slot1,metrics_am_con3{1,2}(slot1),'<','Color',all_colors(2, :));
hold on
con_svd = plot(slot3,metrics_svd_con3{1,2}(slot3),'o','Color',all_colors(3, :));
hold on
con_sota = plot(slot2,metrics_sota_con3{1,2}(slot2),'+','Color',all_colors(5, :));
hold on
con_sbl = plot(slot2,metrics_csbl_con3{1,2}(slot2),'>','Color',all_colors(1, :));
grid on

xlim([0,4000])
ylim([0.00001,10])
yc = yline(metrics_csbl_con3{1,2}(end),'--','NMSE = 0.00054','LineWidth',3);
hold on
ys = yline(metrics_svd_con3{1,2}(end),'--','NMSE = 0.000039','LineWidth',3);
hold on
yso = yline(metrics_sota_con3{1,2}(end),'--','NMSE = 0.00036','LineWidth',3);
hold on
ya = yline(metrics_am_con3{1,2}(end),'--','NMSE = 0.000037','LineWidth',3);
hold on
ys.LabelVerticalAlignment = 'bottom';
yso.LabelVerticalAlignment = 'bottom';
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
title('Known noise variance')
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