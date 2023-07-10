% merge different runs
% two parallel runs, each 50, merge them here
clear
clc
x = load('./Results/compare222_50.mat');
y = load('./Results/compare333_50.mat');

vrs = fieldnames(x);
% Concatenate data
%%
for i = 1:7
 x.('error'){i,2} = cat(3,x.('error'){i,2},y.('error'){i,2});
 x.('time'){i,2} = cat(3,x.('time'){i,2},y.('time'){i,2});
 x.('ser'){i,2} = cat(3,x.('ser'){i,2},y.('ser'){i,2});
end

save('result_paper','-struct','x')