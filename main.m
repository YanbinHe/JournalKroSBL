%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Clear previous work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;

%%% kroSBL code

%%% Import default configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Run the configuration script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('Configuration.m')

func_ctrl.noisy_compare = true;
func_ctrl.noisy_convergence = false;
func_ctrl.noiseless_rank = true;

plt_ctrl.nmse = true;
plt_ctrl.srr = true;

%% Funtional part
%% 1: with the presence of noise, compare the performance of different
% algorithms
if func_ctrl.noisy_compare
    noisy_compare
end
%% 2: with the presence of noise, demonstrate the convergence property
func_ctrl.noisy_convergence = true;
if func_ctrl.noisy_convergence
    noisy_convergence
end
%%
save simu_results_con.mat
%% 3: with the absence of noise, demonstrate the rank refinement
if func_ctrl.noiseless_rank
    noiseless_rank
end
%%
save simu_results.mat
