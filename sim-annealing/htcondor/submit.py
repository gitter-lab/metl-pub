
import uuid
import os
import numpy as np
import pandas as pd
from analysis.preprocessing_gfp_run import preprocess_train_variants_for_sim_anneal,preprocess_only_allow_training
inputs= ['seq2fitness_file',
         'uuid',
         'wild_type',
         'sampler',
         'steps',
         'T0',
         'Tstop',
         'temperature_schedule',
         'number_mutations',
         'mutation_rate',
         'seed',
         'starting_sequence',
         'parallel',
         'AA_constraints',
         'results_dir',
         'split_char',
         'save_final']

defaults = ['seq2fitness',
            str(uuid.uuid4()),
            'skgeelftgvvpilveldgdvnghkfsvsgegegdatygkltlkficttgklpvpwptlvttlsygvqcfsrypdhmkqhdffksampegyvqertiffkddgnyktraevkfegdtlvnrielkgidfkedgnilghkleynynshnvyimadkqkngikvnfkirhniedgsvqladhyqqntpigdgpvllpdnhylstqsalskdpnekrdhmvllefvtaagithgmdelyk',
            'sim_anneal',
            10000,
            1,
            -5,
            'log',
            5,
            3,
            124,
            None,
            1,
            None,
            'results',
            ',',
            'False'
            ]


def generate_args(savepath,**kwargs):
    '''
    for each keyword argument you'd like to modify include a list for that argument to generate files of that type
    output goes to savepath/uuid.txt

    User must make the directory args.

    if kwargs is a list then all lists must be the same length.
    :param kwargs:
    :return: a generated file saved into args with specified run_name as prefix
    '''
    print('stop')

    # do some input checking to see if the
    # input parameters all match in length for different lengths of
    # the lists
    L=[]
    K=[]
    for key,value in zip(kwargs.keys(),kwargs.values()):
        if key not in inputs:
            raise KeyError(f'{key} is not a valid key, accepted keys are: '
                           f'{inputs}')
        if isinstance(value, (list,np.ndarray,tuple,pd.Series)):
            L.append(len(value))
            K.append(key)
        elif isinstance(value,(str,int,float)) or value==None:
            pass
        else:
            raise TypeError(f"unsupported type error for {type(value)} for key : {key},{value},"
                            f"supported types are [str,int,float,list,np.ndarray,tuple,pd.Series,None]")

    # set n_iter for parameter parsing
    n_iter=1
    if len(L)>0:
        assert np.all(L[0]==np.array(L)) , 'the arrays of input with length great '
        n_iter=L[0]


    # initilize the reference
    # dictionary for all the inputs
    args_dict = { }
    for input,default in zip(inputs,defaults):
        if input in kwargs.keys():
            args_dict[input] =kwargs[input]
        else:
            args_dict[input] = default




    # loop over all possible inputs accounting for differences
    # between different input
    for i in range(n_iter):
        # write to the output file a single, have
        # outputs present as a single string
        d = ''
        for key, value in zip(args_dict.keys(),args_dict.values()):
            if key in K:
                d += f"--{key}\n{value[i]}\n"
                if key =='uuid':
                    uuid=value[i]
            elif value==None:
                pass
            else:
                d += f"--{key}\n{value}\n"
                if key=='uuid':
                    uuid= value
        # use uuid as the name of the written file
        with open(os.path.join(savepath,f"{uuid}.txt"),'w') as f:
            f.write(d)



if __name__ == '__main__':
    # AA_constraints = preprocess_train_variants_for_sim_anneal()

    AA_constraints = preprocess_only_allow_training()
    generate_args(os.path.join('args','3d_only_train_5_mutant_10k_run'),
                  uuid= np.arange(10000),
                  seed=np.arange(10000),
                  seq2fitness_file = 'seq2fitness_3d',
                  Tstop=-2,
                  number_mutations =5,
                  mutation_rate=2,
                  steps=10000,
                  AA_constraints= AA_constraints,
                  save_final='True',
                  results_dir='.')

