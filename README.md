# in_silico_opt
This repo is still a work in progress. Feel free to take scripts from here as you need them until we get a finalized version.

Below is a description of each of the files.

- design_tools.py - contains majority of simulated annealing code
- chtc_sa_optimize.sub - chtc submit file to run chtc_sa_optimize.sh for each config in the configs/configs_names.txt file
- chtc_sa_optimize.sh - 
- run_sa.py - python script to run simulated annealing from a config file
- seq2fitness_tools_cnn_0.py - an example of how to set up the seq2fitness file with the init function (called before SA is run in run_sa.py) and the seq2fitness function that given a sequence returns a number
- configs - contains an example for how to set up a config file used in run_sa.py
- output - a directory where your output files should be sent to (this is setup in your config file)

