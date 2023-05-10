import argparse
import importlib
import os
import numpy as np
from utils import AAs
from sampler import *

sample_dict = {'sim_anneal': SimulatedAnnealing, 'random':RandomSampler,'nested':NestedSampler}

def main():

    parser = argparse.ArgumentParser(
        description='python src/run.py seq2fitness.py AbrwbV \n ',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars="@")


    parser.add_argument('--seq2fitness_file', type=str,
                        help='look to example file for reference, must be a function with '
                             'the same name as the filename which returns the fitness'
                             'given a specific sequence')

    parser.add_argument('--uuid',type=str,
                        help='unique identifier as a string for this run')

    parser.add_argument('--wild_type', type=str,
                        help= 'Wildtype argument, all constant wildtype letters should be upper case,'
                              'all mutatable wildtype letters should be lower case i.e. (AymkewG only residues '
                              'at A and G will be held constant in the sampling run)')

    parser.add_argument('--sampler', type=str,
                        help='sampler to run of options simulated annealing (sim_anneal),'
                             ' random walk (random), and nested sampling  (nested),'
                             'default is sim_anneal',
                        default='sim_anneal')

    parser.add_argument('--steps',type=int,
                        help='number of steps to take for the sampler, defaults is 1000',
                        default=1000)


    parser.add_argument('--T0',type=float,
                        help='initial temperature choice, if ',default=3)


    parser.add_argument('--Tstop',type=float,
                         help='the stopping temperature',default=-5)

        
    # need to figure out a way to implement schedule stuff here.

    parser.add_argument('--temperature_schedule',type=str,
                        help='either log or linear',default='log')

    parser.add_argument('--number_mutations',type=int,
                        help='number of allowed mutations to make',default=5)

    parser.add_argument('--mutation_rate',type=int,
                        help='simulated annealing only,'
                             'choose the number of mutations to make to a sequence based off of '
                             'np.random.poisson(mutation_rate) distribution. '
                             'bound the number of mutations to [1,number_mutations-1]',default=3)

    parser.add_argument('--seed',type=int,
                        help='seed for the simulated annealing pseudo-random number generator',default=0)

    parser.add_argument('--starting_sequence',type=str,
                        help='starting mutant, must be same length as wild type and must not contain '
                             'more mutations than the allowed number_mutations,'
                             '--> random is to pick a random mutant',default='random')

    parser.add_argument('--parallel',type=int,
                        help='how many times we should run simulated annealing on this cpu, '
                             'note that this then evaluates the sequences on the fitness function '
                             'on batch mode',default=1)

    parser.add_argument('--AA_constraints', nargs='+',
                        help='--AA_constraints 3GH 4CT '
                             '\n excludes GH at residue 3 (zero indexing) and '
                             '\n CT at residue 4 (zero indexing). '
                             'If you would like to leave an amino acid constant (either its wild type or another)'
                             'at a residue look to the wild_type argument and use an upper case letter. ',
                        type=list,default=[])

    parser.add_argument('--save_final',type=str,
                        help='--save_final True just to save the final output as opposed to the sequences '
                             'from all of the steps.', default='False')

    parser.add_argument('--results_dir',
                        help='results directory',
                        type=str,default='results')

    parser.add_argument('--split_char',help='--split_char between mutant sequences i.e. "T67:U5Y", '
                                            'default is ","',type=str, default=",")


    parsed_args = parser.parse_args()



    print(f'======starting {parsed_args.sampler} ======')
    for key, value in parsed_args.__dict__.items():
        print(key,":", value)
    print('')

    ############### refactor parameters  ######################
    # first define the AA_options and the WT sequence as capital letters
    AA_options = []
    WT = ''
    AA_contraints = parsed_args.AA_constraints
    AA_constraint_dict = {}
    for element in AA_contraints:
        assert len(element) > 1, 'check AA_constraints'
        AA_constraint_dict[int(element[0])] = "".join(element[1:])

    for i, wt_aa in enumerate(parsed_args.wild_type):
        if wt_aa.islower():
            allowed_AAs_at_residue = []
            for AA in AAs:
                if AA_constraint_dict.get(i) != None:
                    if AA not in AA_constraint_dict[i]:
                        allowed_AAs_at_residue.append(AA)
                else:
                    allowed_AAs_at_residue.append(AA)
        else:
            allowed_AAs_at_residue = [wt_aa]

        AA_options.append(allowed_AAs_at_residue)
        WT += wt_aa.upper()

    ########## import seq2fitness library #################
    module_name = parsed_args.seq2fitness_file
    function_name = parsed_args.seq2fitness_file

    try:
        # Attempt to import the specified module
        module = importlib.import_module(module_name)

        # Check if the specified function is in the module
        if not hasattr(module, function_name):
            raise AttributeError(f"The function '{function_name}' is not in the '{module_name}' module.")

    except ImportError:
        # Raise an error if the module is not found
        raise ImportError(f"The module '{module_name}' could not be found.")

    # Get the specified function from the module
    seq2fitness = getattr(module, function_name)



    ##### initilize the sampler to run ######
    sampler= parsed_args.sampler

    assert sampler in sample_dict.keys() , 'invalid sampler input'
    kwargs={}
    if sampler == 'sim_anneal':
        temperature_schedule = parsed_args.temperature_schedule
        assert temperature_schedule in ['log', 'linear'] , 'suport only available for log or linear schedules'
        if temperature_schedule  == 'log':
            handle = np.logspace
        else:
            handle = np.linspace
        mutation_schedule =handle(parsed_args.T0,parsed_args.Tstop,parsed_args.steps)

        #include keyword arguments here specific to just this sampler
        kwargs=  {'number_mutations':parsed_args.number_mutations,'mutation_rate':parsed_args.mutation_rate}

    elif sampler == 'random':
        print('no support currently for a changing mutation schedule')
        mutation_schedule = np.ones((int(parsed_args.steps),))*parsed_args.number_mutations
    elif sampler == 'nested':
        raise Exception('no support for nested sampling yet, sorry!')
    else:
        raise Exception('invalid "--sampler" argument see run.py -h for more details')

    Sampler = sample_dict[parsed_args.sampler](seq2fitness= seq2fitness, AA_options=AA_options,WT= WT,
                                       seed= parsed_args.seed, split_char= parsed_args.split_char,**kwargs)

    # todo : add support for start mutant
    # todo: no support for batch evaluation yet
    df= Sampler.walk(mutation_schedule=mutation_schedule)

    if eval(parsed_args.save_final):
        df=df.tail(1)
    df.to_csv(os.path.join(parsed_args.results_dir,f'{parsed_args.uuid}.csv'),index=False)












if __name__ == '__main__':
    # argument_getter()
    main()