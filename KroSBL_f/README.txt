---------------------------------------------------
code for KroSBL
---------------------------------------------------

# Important information:
Make sure you install tensor_toolbox to be able to start the simulation (already compressed in the directory ./functions, unzip and use).
Put tensor_toolbox under ./functions

# How to start a simulation:

Run main.m

# How to reproduce the figures in the paper
1. unzip results_to_reproduce_figures.zip under ./results
2. 
  a) go to https://drive.google.com/drive/folders/1RpcNHDhvz8oZFQ4FlM63SnJ_eAxToL0N?usp=sharing
  b) download simu_results_con.mat and simu_results.mat, the former is to reproduce the convergence plot, while the latter is to reproduce the rank figure

3. put simu_results_con.mat and simu_results.mat under ./KroSBL_f, run figure_ploting.m, then you should be able to see all the figures

# ./unknown_time_compare is used to reproduce the results in the response letter, you can ignore this part.
