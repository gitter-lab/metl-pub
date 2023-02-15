import numpy as np
from abc import ABC,abstractmethod
import pandas as pd
from tqdm import tqdm



class BaseSampler(ABC):
    def __init__(self,seq2fitness,AA_options,WT,seed,split_char=":"):
        '''
        BaseSampler class for all other methods
        :param seq2fitness: function handle takes in mutations and returns fitness values
        :param AA_options: allowed amino acids at each position
        :param WT: wildtype sequence
        :param seed: seed for random number generator
        :param split_char: split character in mutation data, default is ":"
        '''
        self.seq2fitness=seq2fitness
        self.AA_options=AA_options
        self.WT=WT
        self.rng= np.random.default_rng(seed=seed)
        self.split_char=split_char


    @abstractmethod
    def walk(self, mutation_schedule:list, start_mut=None):
        '''
        parent function of all the sampling methods
        :param mutation_schedule: list of a temperature gradient (simulated annealing) or number
        of mutations to make at each step (random walk), or the percentage of mutations to accept
        between each monte carlo run. (nested sampling)
        :param start_mut: if want to start mutations at a specific spot
        :return: a pandas dataframe containing all relevant information from the run
        '''
        pass


    def generate_all_point_mutants(self):
        """Generate all possible single point mutants of a sequence from AA_options"""
        all_mutants = []
        for pos in range(len(self.WT)):
            for aa in self.AA_options[pos]:
                if self.WT[pos] != aa:  # not mutating to self
                    mut = self.WT[pos] + str(pos) + aa
                    all_mutants.append(mut)
        return all_mutants

    def mut2seq(self, mutations):
        '''
        Create mutations in form of A94T to a mutant sequence from wild type
        :param mutations: list of mutations in form of ["A94T", "H99R"] or "A94T <self.split_char> H99R"
        :return: mutated sequence from wild type
        '''

        if type(mutations) is str:
            mutations = mutations.split(self.split_char)
        mutant_seq=list(self.WT)
        for mut in mutations:
            wt_aa,pos,new_aa =mut[0],int(mut[1:-1]),mut[-1]
            assert wt_aa == self.WT[pos], 'Warning: WT residue in mutation %s does not match WT sequence'
            mutant_seq[pos] = new_aa
        return ''.join(mutant_seq)
    def generate_mutant(self, all_mutants:list, number_of_mutations:int, current_mutant=None):
        '''
        generate mutant from current mutant
        :param all_mutants: all mutant library of all allowed mutations.
        This can contain mutations that are in current mutant. However, it cannot contain mutations that are never allowed,
        such as those specified by AA_options
        :param number_of_mutations: number of mutations to allow to current mutant, if current mutant is specified
        it will append only until that amount of mutations.
        :param current_mutant: starting mutant , default is no mutant, and it will generate a random
        mutant with specified number_of_mutations
        :return: generated mutant of the form T5I,U5L where the ',' is the specified split character

        Example:
        new_mutant= generate_mutant(all_mutants=['S0A', 'S0C', 'S0D', 'S0E', 'S0F', 'S0G', 'S0H',...],
                                    number_of_mutations=5,
                                    current_mutant=['S0E','I6Y','L9P'])
        print(new_mutant)
        **** display:  'S0E,I6Y,L9P,N18T,K1A'
        '''

        if current_mutant is None:
            current_mutant = []

        assert number_of_mutations >= len(current_mutant), \
                'number of mutations must be greater than current mutant'

        # find the positions which are not allowed
        not_allowed_positions = [int(mutant[1:-1]) for mutant in current_mutant]
        generated_mutant = list(current_mutant)
        for _ in range(number_of_mutations-len(current_mutant)):
            updated_mutations = [mutant for mutant in all_mutants if
                                 int(mutant[1:-1]) not in not_allowed_positions]
            mutation_choice = self.rng.choice(updated_mutations, size=1)[0]
            not_allowed_positions.append(int(mutation_choice[1:-1]))
            generated_mutant.append(mutation_choice)

        assert np.unique([int(mutant[1:-1]) for mutant in generated_mutant]).shape[0] == number_of_mutations,\
                            'the number of unique positions in current mutant must equal that of the suggested number ' \
                            'of mutations'
        assert np.all([mutant[-1] != self.WT[int(mutant[1:-1])] for mutant in generated_mutant]),\
                'a suggested mutation mutation is going to wildtype, try checking all_mutants'

        return f'{self.split_char}'.join(generated_mutant)

class RandomSampler(BaseSampler):
    ''''
    RandomSampler that generates the number of suggested random mutations
    as given by the mutation schedule, all mutations are taken from wild type.
    '''
    def walk(self, mutation_schedule, start_mut=None, ):
        # todo: how to handle start mutation here
        all_mutants=self.generate_all_point_mutants()
        M,F =[],[]
        for num_mut in mutation_schedule:
            mutant= self.generate_mutant(all_mutants,num_mut)
            fitness = self.seq2fitness(mutant)
            F.append(fitness)

        return pd.DataFrame(data=[M,F],columns=['mutant','fitness'])

class SimulatedAnnealing(BaseSampler):
    def __init__(self, seq2fitness, AA_options, WT, seed,split_char,**kwargs):
        super().__init__(seq2fitness, AA_options, WT, seed,split_char)

        self.mutation_rate=kwargs['mutation_rate']
        self.number_mutations=kwargs['number_mutations']

    def walk(self, mutation_schedule, start_mut=None):
        all_mutation_library = self.generate_all_point_mutants()
        if start_mut is None:
            start_mut = self.generate_mutant(all_mutation_library,self.number_mutations)

        fit = self.seq2fitness(start_mut,self.WT)

        # choose the number of mutations to make to current sequence
        # BM = Best Mutant  , BF= Best Fitness
        # M =current mutant,  F = current fitness
        # n = number of mutations to make for a given step
        # T = temperature schedule (which is inherently the number of steps

        N,M,F,BM,BF= [len(start_mut.split(self.split_char))],[start_mut],[fit],[start_mut],[fit]

        for T in tqdm(mutation_schedule):

            current_mutant = M[-1]
            current_fitness = F[-1]

            n = self.rng.poisson(self.mutation_rate)
            n = min([self.number_mutations - 1, max([1, n])])  # bound n within the range [1,number_mututations-1]
            N.append(n)
            reduced_mutant= self.rng.choice(current_mutant.split(self.split_char),
                                            size=self.number_mutations-n,replace=False)

            temp_mutant = self.generate_mutant(all_mutants=all_mutation_library,
                                               number_of_mutations=self.number_mutations,
                                               current_mutant=reduced_mutant)

            assert len(temp_mutant.split(self.split_char))==len(current_mutant.split(self.split_char)),\
                'length of current mutant and temp mutant must be equal'

            temp_fitness= self.seq2fitness(temp_mutant,self.WT)

            # now check to see if it is the best mutant
            if temp_fitness > BF[-1]:
                BF.append(temp_fitness)
                BM.append(temp_mutant)
            else:
                BF.append(BF[-1])
                BM.append(BM[-1])


            # Simulated annealing acceptance criteria:
            #   If mutant is better than current seq, move to mutant seq
            #   If mutant is worse than current seq, move to mutant with some
            #   exponentially decreasing probability with delta_F
            delta_F = temp_fitness - F[-1]  # new seq is worse if delta_F is neg.
            # print(np.exp(min([0,delta_F/(10*T)])))
            if np.exp(min([0, delta_F / (T)])) > self.rng.random():
                M.append(temp_mutant)
                F.append(temp_fitness)
            else:
                M.append(current_mutant)
                F.append(current_fitness)

        # start mutant didn't have any temperature
        mutation_schedule=np.hstack(([None],mutation_schedule))
        df= pd.DataFrame(data=[M, F, BM, BF, N, mutation_schedule]).T
        df.columns=['current_mutant','current_fitness','best_mutant','best_fitness',
                                     'number_of_mutations','temperature']

        # return all the data found from the run
        return df




class BatchSimulatedAnnealing(BaseSampler):
    def __init__(self, seq2fitness, AA_options, WT,seed):
        super().__init__(seq2fitness, AA_options, WT,seed)



# class HillClimber(BaseSampler):
# '''
# Greedy HillClimber
# '''


class NestedSampler(BaseSampler):
    def __init__(self, seq2fitness, AA_options, WT, seed):
        super().__init__(seq2fitness, AA_options, WT, seed)


       # get scores for all mutants



