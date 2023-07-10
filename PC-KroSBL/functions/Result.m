AVG = 20;
SNRl = [5,10,15,20,25,30];
overheadi = [4,12];
%% ii-overhead jj-snr
load re1.mat
II = 2;
JJ = 6;
results1_kro1 = zeros(II,JJ);
results1_kro2 = zeros(II,JJ);
results1_cla = zeros(II,JJ);
results1_omp = zeros(II,JJ);
results1_r11 = zeros(II,JJ);
results1_r12 = zeros(II,JJ);
results1_r2 = zeros(II,JJ);
results1_ro = zeros(II,JJ);
results1_s11 = zeros(II,JJ);
results1_s12 = zeros(II,JJ);
results1_s2 = zeros(II,JJ);
results1_sp = zeros(II,JJ);
results1_so = zeros(II,JJ);

for ii = 1:II % for each overhead situation
    for jj = 3:8 % snr
        % error
        results1_kro1(ii,jj-2) = sum(errorl1(ii,jj,:))/AVG;
        results1_kro2(ii,jj-2) = sum(errorl2(ii,jj,:))/AVG;
        results1_cla(ii,jj-2) = sum(errorc(ii,jj,:))/AVG;
        results1_omp(ii,jj-2) = sum(er(ii,jj,:))/AVG;
        % recovery rate
        results1_r11(ii,jj-2) = sum(r11(ii,jj,:))/AVG;
        results1_r12(ii,jj-2) = sum(r12(ii,jj,:))/AVG;
        results1_r2(ii,jj-2) = sum(r2(ii,jj,:))/AVG;
        results1_ro(ii,jj-2) = sum(ro(ii,jj,:))/AVG;
        % ser
        results1_s11(ii,jj-2) = sum(s11(ii,jj,:))/AVG;
        results1_s12(ii,jj-2) = sum(s12(ii,jj,:))/AVG;
        results1_s2(ii,jj-2) = sum(s2(ii,jj,:))/AVG;
        results1_sp(ii,jj-2) = sum(sperf(ii,jj,:))/AVG;
        results1_so(ii,jj-2) = sum(so(ii,jj,:))/AVG;
    end
end
%% ii-overhead jj-snr
load re2.mat
II = 2;
JJ = 6;
results2_kro1 = zeros(II,JJ);
results2_kro2 = zeros(II,JJ);
results2_cla = zeros(II,JJ);
results2_omp = zeros(II,JJ);
results2_r11 = zeros(II,JJ);
results2_r12 = zeros(II,JJ);
results2_r2 = zeros(II,JJ);
results2_ro = zeros(II,JJ);
results2_s11 = zeros(II,JJ);
results2_s12 = zeros(II,JJ);
results2_s2 = zeros(II,JJ);
results2_sp = zeros(II,JJ);
results2_so = zeros(II,JJ);

for ii = 1:II % for each overhead situation
    for jj = 1:JJ % snr
        % error
        results2_kro1(ii,jj) = sum(errorl1(ii,jj,:))/AVG;
        results2_kro2(ii,jj) = sum(errorl2(ii,jj,:))/AVG;
        results2_cla(ii,jj) = sum(errorc(ii,jj,:))/AVG;
        results2_omp(ii,jj) = sum(er(ii,jj,:))/AVG;
        % recovery rate
        results2_r11(ii,jj) = sum(r11(ii,jj,:))/AVG;
        results2_r12(ii,jj) = sum(r12(ii,jj,:))/AVG;
        results2_r2(ii,jj) = sum(r2(ii,jj,:))/AVG;
        results2_ro(ii,jj) = sum(ro(ii,jj,:))/AVG;
        % ser
        results2_s11(ii,jj) = sum(s11(ii,jj,:))/AVG;
        results2_s12(ii,jj) = sum(s12(ii,jj,:))/AVG;
        results2_s2(ii,jj) = sum(s2(ii,jj,:))/AVG;
        results2_sp(ii,jj) = sum(sperf(ii,jj,:))/AVG;
        results2_so(ii,jj) = sum(so(ii,jj,:))/AVG;
    end
end
%% ii-overhead jj-snr
load re3.mat
II = 2;
JJ = 6;
results3_kro1 = zeros(II,JJ);
results3_kro2 = zeros(II,JJ);
results3_cla = zeros(II,JJ);
results3_omp = zeros(II,JJ);
results3_r11 = zeros(II,JJ);
results3_r12 = zeros(II,JJ);
results3_r2 = zeros(II,JJ);
results3_ro = zeros(II,JJ);
results3_s11 = zeros(II,JJ);
results3_s12 = zeros(II,JJ);
results3_s2 = zeros(II,JJ);
results3_sp = zeros(II,JJ);
results3_so = zeros(II,JJ);

for ii = 1:II % for each overhead situation
    for jj = 1:JJ % snr
        % error
        results3_kro1(ii,jj) = sum(errorl1(ii,jj,:))/AVG;
        results3_kro2(ii,jj) = sum(errorl2(ii,jj,:))/AVG;
        results3_cla(ii,jj) = sum(errorc(ii,jj,:))/AVG;
        results3_omp(ii,jj) = sum(er(ii,jj,:))/AVG;
        % recovery rate
        results3_r11(ii,jj) = sum(r11(ii,jj,:))/AVG;
        results3_r12(ii,jj) = sum(r12(ii,jj,:))/AVG;
        results3_r2(ii,jj) = sum(r2(ii,jj,:))/AVG;
        results3_ro(ii,jj) = sum(ro(ii,jj,:))/AVG;
        % ser
        results3_s11(ii,jj) = sum(s11(ii,jj,:))/AVG;
        results3_s12(ii,jj) = sum(s12(ii,jj,:))/AVG;
        results3_s2(ii,jj) = sum(s2(ii,jj,:))/AVG;
        results3_sp(ii,jj) = sum(sperf(ii,jj,:))/AVG;
        results3_so(ii,jj) = sum(so(ii,jj,:))/AVG;
    end
end
%% ii-overhead jj-snr
num = 3;
for ii = 1:II % for each overhead situation
    for jj = 1:JJ % snr
        % error
        results_kro1(ii,jj) = (results1_kro1(ii,jj)+results2_kro1(ii,jj)+results3_kro1(ii,jj))/num;
        results_kro2(ii,jj) = (results1_kro2(ii,jj)+results2_kro2(ii,jj)+results3_kro2(ii,jj))/num;
        results_cla(ii,jj) = (results1_cla(ii,jj)+results2_cla(ii,jj)+results3_cla(ii,jj))/num;
        results_omp(ii,jj) = (results1_omp(ii,jj)+results2_omp(ii,jj)+results3_omp(ii,jj))/num;
        % recovery rate
        results_r11(ii,jj) = (results1_r11(ii,jj)+results2_r11(ii,jj)+results3_r11(ii,jj))/num;
        results_r12(ii,jj) = (results1_r12(ii,jj)+results2_r12(ii,jj)+results3_r12(ii,jj))/num;
        results_r2(ii,jj) = (results1_r2(ii,jj)+results2_r2(ii,jj)+results3_r2(ii,jj))/num;
        results_ro(ii,jj) = (results1_ro(ii,jj)+results2_ro(ii,jj)+results3_ro(ii,jj))/num;
        % ser
        results_s11(ii,jj) = (results1_s11(ii,jj)+results2_s11(ii,jj)+results3_s11(ii,jj))/num;
        results_s12(ii,jj) = (results1_s12(ii,jj)+results2_s12(ii,jj)+results3_s12(ii,jj))/num;
        results_s2(ii,jj) = (results1_s2(ii,jj)+results2_s2(ii,jj)+results3_s2(ii,jj))/num;
        results_sp(ii,jj) = (results1_sp(ii,jj)+results2_sp(ii,jj)+results3_sp(ii,jj))/num;
        results_so(ii,jj) = (results1_so(ii,jj)+results2_so(ii,jj)+results3_so(ii,jj))/num;
    end
end
%%
% figure
% for jj = 1:6 % for different SNR situation
% plot(overheadi,results_kro1(:,jj),'-ro','DisplayName','KroSBL');
% hold on
% plot(overheadi,results_kro2(:,jj),'-bx','DisplayName','SBL');
% hold on
% plot(overheadi,results_cla(:,jj),'-b^','DisplayName','SBL');
% hold on
% plot(overheadi,results_omp(:,jj),'-kv','DisplayName','OMP');
% hold on
% end
% xlabel('#Overhead')
% ylabel('NMSE')
% grid on
% set(gca, 'yscale', 'log');
% set(0,'DefaultLineLineWidth',3)
% set(0,'DefaultAxesFontSize',16)
% set(0,'DefaultLineMarkerSize',14)
% set(0,'DefaultAxesFontWeight','bold')
% set(gca,'FontSize',16)
% set(get(gca,'Xlabel'),'FontSize',16)
% set(get(gca,'Ylabel'),'FontSize',16)
% set(get(gca,'Title'),'FontSize',16)
% set(get(gca,'Xlabel'),'FontWeight','bold')
% set(get(gca,'Ylabel'),'FontWeight','bold')
% set(get(gca,'Title'),'FontWeight','bold')
%%
figure
plot(SNRl,results_kro1(1,:),'-.o','Color',[252 170 103]/255,'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_kro2(1,:),'-.x','Color',[189 30 30]/255,'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_cla(1,:),'-.^','Color',[0 70 222]/255,'DisplayName','SBL');
hold on
plot(SNRl,results_omp(1,:),'-.v','Color',[71 51 53]/255,'DisplayName','OMP');
hold on

plot(SNRl,results_kro1(2,:),'-o','Color',[252 170 103]/255,'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_kro2(2,:),'-x','Color',[189 30 30]/255,'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_cla(2,:),'-^','Color',[0 70 222]/255,'DisplayName','SBL');
hold on
plot(SNRl,results_omp(2,:),'-v','Color',[71 51 53]/255,'DisplayName','OMP');
hold on

xlabel('SNR')
ylabel('NMSE')
grid on
set(gca, 'yscale', 'log');
set(gca, 'yscale', 'log');
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',16)
set(get(gca,'Xlabel'),'FontSize',16)
set(get(gca,'Ylabel'),'FontSize',16)
set(get(gca,'Title'),'FontSize',16)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
%%
figure
for ii = 1 % for each SNR case
plot(SNRl,results_r11(ii,:),'-.o','Color',[252 170 103]/255,'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_r12(ii,:),'-.x','Color',[189 30 30]/255,'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_r2(ii,:),'-.^','Color',[0 70 222]/255,'DisplayName','SBL');
hold on
end
for ii = 2 % for each SNR case
plot(SNRl,results_r11(ii,:),'-o','Color',[252 170 103]/255,'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_r12(ii,:),'-x','Color',[189 30 30]/255,'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_r2(ii,:),'-^','Color',[0 70 222]/255,'DisplayName','SBL');
hold on
end
xlabel('SNR')
ylabel('Recover rate')
grid on
% set(gca, 'yscale', 'log');
% set(gca, 'yscale', 'log');
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',16)
set(get(gca,'Xlabel'),'FontSize',16)
set(get(gca,'Ylabel'),'FontSize',16)
set(get(gca,'Title'),'FontSize',16)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
%%
figure
for ii = 1
plot(SNRl,results_s11(ii,:),'-.o','Color',[252 170 103]/255,'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_s12(ii,:),'-.x','Color',[189 30 30]/255,'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_s2(ii,:),'-.^','Color',[0 70 222]/255,'DisplayName','SBL');
hold on
plot(SNRl,results_sp(ii,:),'-.d','Color',[62 43 109]/255,'DisplayName','Perfect CSI');
hold on
plot(SNRl,results_so(ii,:),'-.v','Color',[71 51 53]/255,'DisplayName','OMP');
hold on
end

for ii = 2
plot(SNRl,results_s11(ii,:),'-o','Color',[252 170 103]/255,'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_s12(ii,:),'-x','Color',[189 30 30]/255,'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_sp(ii,:),'-^','Color',[0 70 222]/255,'DisplayName','SBL');
hold on
plot(SNRl,results_sp(ii,:),'-d','Color',[62 43 109]/255,'DisplayName','Perfect CSI');
hold on
plot(SNRl,results_so(ii,:),'-v','Color',[71 51 53]/255,'DisplayName','OMP');
hold on
end

xlabel('SNR')
ylabel('SE')
grid on
%%
% figure
% for ii = 1:1
% plot(SNRl,reshape(time1(1,ii,:)/AVG,JJ,1),'-.o','Color',[252 170 103]/255,'DisplayName','KroSBL with SVD');
% hold on
% plot(SNRl,reshape(time2(1,ii,:)/AVG,JJ,1),'-.x','Color',[189 30 30]/255,'DisplayName','KroSBL with Alternating');
% hold on
% plot(SNRl,reshape(time(2,ii,:)/AVG,JJ,1),'-.^','Color',[0 70 222]/255,'DisplayName','SBL');
% hold on
% plot(SNRl,reshape(time(3,ii,:)/AVG,JJ,1),'-.v','Color',[71 51 53]/255,'DisplayName','OMP');
% hold on
% end
% set(gca, 'yscale', 'log');
% 
% for ii = 2:2
% plot(SNRl,reshape(time1(1,ii,:)/AVG,JJ,1),'-o','Color',[252 170 103]/255,'DisplayName','KroSBL with SVD');
% hold on
% plot(SNRl,reshape(time2(1,ii,:)/AVG,JJ,1),'-x','Color',[189 30 30]/255,'DisplayName','KroSBL with Alternating');
% hold on
% plot(SNRl,reshape(time(2,ii,:)/AVG,JJ,1),'-^','Color',[0 70 222]/255,'DisplayName','SBL');
% hold on
% plot(SNRl,reshape(time(3,ii,:)/AVG,JJ,1),'-v','Color',[71 51 53]/255,'DisplayName','OMP');
% hold on
% end
% set(gca, 'yscale', 'log');
% xlabel('SNR')
% ylabel('Seconds')
% grid on