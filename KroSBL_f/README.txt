---------------------------------------------------
Simulator for Parallel-Concatenated Turbo Codes
---------------------------------------------------
(c) Dr. Christoph Studer 2011 (studer@rice.edu)
---------------------------------------------------

# Important information:

If you are thinking of contacting us, please do not e-mail the author to ask for download instructions, installation guidelines, or the toolbox itself. Note that we will NOT help to debug user-generated code that is not included in the provided package. If, however, you notice a bug in our code, please be so kind to contact the author of this software: C. Studer (studer@rice.edu). 

The package is supplied "as is", without any accompanying support services, maintenance, or future updates. We make no warranties, explicit or implicit, that the software contained in this package is free of error or that it will meet your requirements for any particular application. It should not be relied on for any purpose where incorrect results could result in loss of property, personal injury, liability or whatsoever. Remember: If you do use our software for any such purpose, it is at your own risk. The authors disclaim all liability of any kind, either direct or consequential, resulting from your use of these programs.

# How to start a simulation:

Note that Matlab must be started in the main folder in order to automatically include all paths defined in the userpath.m file. Note that the necessary paths can also be included manually (using the Matlab GUI). Starting a simulation is straightforward: The simulator bases on parameter files (found in the param/ folder), which are used to define all necessary simulation settings and also start the turbo encoding and decoding simulations. For example, type

>> ERR_PCTC_LTE_V3_R12_1952b_I10(0) 

which starts a simulation with a parallel-concatenated turbo code (PCTC) in an AWGN channel using a block-length 1952, code rate 1/2 (punctured), 10 inner decoder iterations, and the 3GPP LTE interleaver. The decoder algorithm can be specified in the parameter file and either performs the sum-product algorithm (SPA) or the max-product algorithm (MPA, also known as max-log approximation). The value 0 determines the random seed and is useful in combination with high-throughput computing, e.g., in combination with Condor. 

The folder "codes" contains a few example codes (e.g., for 3GPP LTE or the code originally used by Berrou et al.). For each PCTC, there exists a Matlab script describing the code properties, e.g., the file PCTC_LTE_V3_R12_1952b.m contains all the essential information about the length-1952 rate-1/2 PCTC code used in the example above. Before a simulation with a new code can be used, one needs to execute this script (which calls genmat.m) in the codes/ folder, which results in a corresponding .mat file, e.g., PCTC_LTE_V3_R12_1952b.mat (also located in the codes/ folder). This .mat file is then used by the parameter files and used by the simulator.

Remember: Since the SPA and MPA decoding algorithms are written in C (to reduce the simulation time significantly), it might be necessary that you have to recompile BCJR.c and BCJRopt.c for your target platform. 

Important: We highly recommend you to execute the code step-by-step (using Matlab's debug mode) in order to gain some understanding of the simulator and how PCTCs are described, constructed, and stored. 

# Version 1.0 (Dec 11, 2011) - initial version for public access