#!/usr/bin/env python
# coding: utf-8

import design_tools as tools
import pickle
import random
import numpy as np
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
    AA_options = [tuple([AA for AA in AAs]) for i in range(len(config['WT']))]
    AA_options.pop(0)
    AA_options.insert(0, ['M'])

    seq2fitness_tools = importlib.__import__(config['seq2fitness_tools'])
    seq2fitness_handler = seq2fitness_tools.seq2fitness_handler()

    sa_optimizer = tools.SA_optimizer(seq2fitness_handler.seq2fitness, config['WT'], AA_options,
            config['num_mut'], mut_rate=config['mut_rate'], nsteps=config['nsteps'],
            cool_sched=config['cool_sched'])
    best_mut, fitness = sa_optimizer.optimize(seed=config['seed'])
    with open(config['export_best_seqs'], 'wb') as f:
        pickle.dump([best_mut, fitness], f)

    if config['save_plot_trajectory']:
        sa_optimizer.plot_trajectory(savefig_name=config['file_plot_trajectory'])


if __name__ == '__main__':
    config = load_config(sys.argv[1])
    run_simulated_annealing(config)

