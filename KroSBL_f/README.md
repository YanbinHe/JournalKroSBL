---------------------------------------------------
code for KroSBL
---------------------------------------------------

# Important information:
1. Make sure you install tensor_toolbox to be able to start the simulation (already compressed in the directory ./functions, unzip and use).
Put tensor_toolbox under ./functions
2. A run for a simulation (200 average) might take 25-30 hours

# How to start a simulation:

Run main.m

# How to reproduce the figures in the paper
1. create a new directory under main ./results
2. go to https://drive.google.com/drive/folders/1RpcNHDhvz8oZFQ4FlM63SnJ_eAxToL0N?usp=sharing
3. download, put, and unzip noise_compare.zip under ./results
4. run figure_ploting.m
5. you will see Fig. 3.
6. download, put, and unzip convergence.zip under ./results under ./KroSBL_f, run figure_ploting.m, then you should be able to see Figs. 2 and 6.
7. download, put, and unzip unknown_nsq_nmse_m.zip under ./unknwon_time_compare/results, run the function at the end of main.m, you will see Fig. 5.
8. run review_convergence.m, you will see Fig. 4.

# All the files are tested in Matlab 2021b, but any later version should also work.
