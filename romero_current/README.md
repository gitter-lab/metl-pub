# in_silico_opt

**Bryce accidentaly deleted a bunch of stuff, look to main branch for most current version. This is just for**
**Bryce's easy reference.**

This repo is still a work in progress. Feel free to take scripts from here as you need them until we get a finalized version.

Below is a description of each of the files.

- Makefile - contains commands to run simulated annealing on chtc. You will need to modify what is included in the make sa_optimize_files.tar.gz command to include the files and directories relevant to your code. Example is below
- design_tools.py - contains majority of simulated annealing code
- chtc_sa_optimize.sub - chtc submit file to run chtc_sa_optimize.sh for each config in the configs/configs_names.txt file
- chtc_sa_optimize.sh - bash script run on chtc that takes config files as passes config file
- run_sa.py - python script to run simulated annealing from a config file
- seq2fitness_tools_cnn_0.py - an example of how to set up the seq2fitness file with the init function (called before SA is run in run_sa.py) and the seq2fitness function that given a sequence returns a number
- configs - contains an example for how to set up a config file used in run_sa.py
- output - a directory where your output files should be sent to (this is setup in your config file)

## Example of how to run/what to change:
- Create your own seq2fitness_tools python script
    - import relevant packages
    - write the init function (initialize any models/variables required for you seq2fitness function, in particular use self.your_var to keep track of relevant variables)
    - write the seq2fitness function (take the sequence that is passed and return a number corresponding to the fitness, use self.your_var to reference any variables that you created in the init function)
- Create an chtc-compatible environment for your seq2fitness tools
    - be sure to add matplotlib, numpy, and yaml in your environment
    - add this environment to the .sub file in place of the docker image
- Create configs for each of your SA runs (can use a python script to do this)
    - add the name of the config file(s) to the configs_names.txt file, each on separate lines
- Update the Makefile to include any files and directories to your .tar.gz file that you will submit to chtc
- Submit your job to CHTC
  - from the in_silico_opt folder run the make command: <code>make</code>
  - This will automatically generate/update the .tar.gz file for submission and submit all of the jobs from the configs/configs_names.txt file


