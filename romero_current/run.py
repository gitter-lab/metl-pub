#!/usr/bin/env python
# coding: utf-8

import design_tools as tools
import pickle
import random
import numpy as np
import pandas as pd
import sys
import yaml
import importlib

AAs = 'ACDEFGHIKLMNPQRSTVWY'

def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def run_simulated_annealing(config):

    # constants maybe a parameter to define the backbone of the sequence.
    # or have one input parameter
    AA_options = [tuple([AA for AA in AAs]) for i in range(len(config['WT']))]
    AA_options.pop(0)
    AA_options.insert(0, ('M')) # this sets the first mutation as a constant 'M' ...

    # define the sequence to fitness handler
    seq2fitness_tools = importlib.__import__(config['seq2fitness_tools'])
    seq2fitness_handler = seq2fitness_tools.seq2fitness_handler()


    sa_optimizer = tools.SA_optimizer(seq2fitness_handler.seq2fitness, config['WT'], AA_options,
            config['num_mut'], mut_rate=config['mut_rate'], nsteps=config['nsteps'],
            cool_sched=config['cool_sched'])

    # return fitness trajectory
    fitness_trajectory = sa_optimizer.optimize(seed=config['seed'])

    # input all the run parameters
    
    if config['save_plot_trajectory']:
        sa_optimizer.plot_trajectory(savefig_name=config['file_plot_trajectory'])


if __name__ == '__main__':
    # config = load_config(sys.argv[1])
    config = load_config('../configs/config_num_mut_5_nstep_50000_sim_200_seq2fitness_tools_cnn_0.txt')
    run_simulated_annealing(config)

