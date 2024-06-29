
import numpy as np
def seq2fitness_example(mutant:str,WT):
    '''
    my example sequence to fitness function
    :param mutant: string of mutants, "," as specified by arguments "E3K,G102S"
    :param WT: Wildtype sequence
    :return: a scalar value of fitness
    '''

    # note your split character must align with that specified in args document
    return np.random.poisson(len(mutant.split(',')))
