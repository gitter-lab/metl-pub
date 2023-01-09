import numpy as np
import random
import matplotlib.pyplot as plt

# TODO: decide if index from 0 or 1 for mutations


def generate_all_point_mutants(seq, AA_options=None):
    """Generate all possible single point mutants of a sequence
    Arguments:
        seq: starting seq
        AA_options: list of amino acid options at each position, if none defaults to all 20 AAs (default None)
    """
    all_mutants = []
    for pos in range(len(seq)):
        for aa in AA_options[pos]:
            if seq[pos] != aa:  # not mutating to self
                mut = seq[pos] + str(pos) + aa
                all_mutants.append(mut)
    return all_mutants


def mut2seq(seq, mutations):
    """Create mutations in form of A94T to seq
    Arguments:
        seq: starting seq
        mutations: list of mutations in form of ["A94T", "H99R"] or "A94T,H99R"
    """
    mutant_seq = seq
    if type(mutations) is str:
        mutations = mutations.split(',')
    for mut in mutations:
        pos = int(mut[1:-1])
        newAA = mut[-1]
        if mut[0] != seq[pos]:
            print('Warning: WT residue in mutation %s does not match WT sequence' % mut)
        mutant_seq = mutant_seq[:pos] + newAA + mutant_seq[pos + 1:]
    return mutant_seq


def find_top_n_mutations(seq2fitness, all_mutants, WT, n=10):
    # evaluate fitness of all single mutants from WT
    single_mut_fitness = []
    for mut in all_mutants:
        seq = mut2seq(WT, (mut,))
        fit = seq2fitness(seq)
        single_mut_fitness.append((mut, fit))

    # find the best mutation per position
    best_mut_per_position = []
    for pos in range(len(WT)):
        best_mut_per_position.append(max([m for m in single_mut_fitness if int(m[0][1:-1]) == pos], key=lambda x: x[1]))

    # take the top n
    sorted_by_fitness = sorted(best_mut_per_position, key=lambda x: x[1], reverse=True)
    topn = [m[0] for m in sorted_by_fitness[:n]]
    topn = tuple([n[1] for n in sorted([(int(m[1:-1]), m) for m in topn])])  # sort by position

    return topn


def generate_random_mut(WT, AA_options, num_mut):
    # Want to make the probability of getting any mutation the same:
    AA_mut_options = []
    for WT_AA, AA_options_pos in zip(WT, AA_options):
        if WT_AA in AA_options_pos:
            options = list(AA_options_pos).copy()
            options.remove(WT_AA)
            AA_mut_options.append(options)
    mutations = []

    for n in range(num_mut):

        num_mut_pos = sum([len(row) for row in AA_mut_options])
        prob_each_pos = [len(row) / num_mut_pos for row in AA_mut_options]
        rand_num = random.random()
        for i, prob_pos in enumerate(prob_each_pos):
            rand_num -= prob_pos
            if rand_num <= 0:
                mutations.append(WT[i] + str(i) + random.choice(AA_mut_options[i]))
                AA_mut_options.pop(i)
                AA_mut_options.insert(i, [])
                break
    return ','.join(mutations)


class SA_optimizer:
    def __init__(self, seq2fitness, WT, AA_options, num_mut, mut_rate=1, nsteps=1000, cool_sched='log'):
        self.seq2fitness = seq2fitness
        self.WT = WT
        self.AA_options = AA_options
        self.num_mut = num_mut
        self.mut_rate = mut_rate
        self.nsteps = nsteps
        self.cool_sched = cool_sched

    def optimize(self, start_mut=None, seed=0):
        random.seed(0)

        if start_mut is None:
            start_mut = generate_random_mut(self.WT, self.AA_options, self.num_mut).split(',')
        all_mutants = generate_all_point_mutants(self.WT, self.AA_options)

        assert ((self.cool_sched == 'log') or (self.cool_sched == 'lin')), 'cool_sched must be \'log\' or \'lin\''
        if self.cool_sched == 'log':
            temp = np.logspace(3, -5, self.nsteps)
        if self.cool_sched == 'lin':
            temp = np.linspace(1000, 1e-9, self.nsteps)

        print('Simulated Annealing Progress: ')
        seq = mut2seq(self.WT, start_mut)
        fit = self.seq2fitness(seq)
        current_seq = [start_mut, fit]
        self.best_seq = [start_mut, fit]
        self.fitness_trajectory = [[fit, fit]]  # fitness_trajectory = [best_fit, current_fit]

        # for loop over decreasing temperatures
        for T in temp:
            mutant = list(current_seq[0])

            # choose the number of mutations to make to current sequence
            n = np.random.poisson(self.mut_rate)
            n = min([self.num_mut - 1, max([1, n])])  # bound n within the range [1,num_mut-1]

            # remove random mutations from current sequence until it contains (num_mut-n) mutations
            while len(mutant) > (self.num_mut - n):
                mutant.pop(random.choice(range(len(mutant))))

            # add back n random mutations to generate new mutant
            occupied = [m[1:-1] for m in mutant]  # postions that already have a mutation
            mut_options = [m for m in all_mutants if m[1:-1] not in occupied]  # mutations at unoccupied positions
            while len(mutant) < self.num_mut:
                mutant.append(random.choice(mut_options))
                occupied = [m[1:-1] for m in mutant]
                mut_options = [m for m in all_mutants if m[1:-1] not in occupied]

            # sort mutations by position to clean up
            mutant = tuple([n[1] for n in sorted([(int(m[1:-1]), m) for m in mutant])])

            # evaluate fitness of new mutant
            fitness = self.seq2fitness(mut2seq(self.WT, mutant))

            # if mutant is better than best sequence, reassign best sequence
            if fitness > self.best_seq[1]:
                self.best_seq = [mutant, fitness]
                # print('fitness = %0.2f, mutations = %s' % (fitness,''.join([i.ljust(5) for i in mutant])))

            # Simulated annealing acceptance criteria:
            #   If mutant is better than current seq, move to mutant seq
            #   If mutant is worse than current seq, move to mutant with some exponentially decreasing probability with delta_F
            delta_F = fitness - current_seq[1]  # new seq is worse if delta_F is neg.
            # print(np.exp(min([0,delta_F/(10*T)])))
            if np.exp(min([0, delta_F / (T)])) > random.random():
                current_seq = [mutant, fitness]

            # store the current fitness
            self.fitness_trajectory.append([self.best_seq[1], current_seq[1]])

        return self.best_seq  # returns [best_mut, best_fit]

    def plot_trajectory(self, savefig_name=None):
        plt.plot(np.array(self.fitness_trajectory)[:, 0])
        plt.plot(np.array(self.fitness_trajectory)[:, 1])
        plt.xlabel('Step')
        plt.ylabel('Fitness')
        plt.legend(['Best mut found', 'Current mut'])
        if savefig_name is None:
            plt.show()
        else:
            plt.savefig(savefig_name)
        plt.close()


class Hill_climber():

    def __init__(self, seq2fitness, WT, AA_options, num_mut, mut_rate=1, num_restarts=100, max_steps=1000,
                 seq2fitness_many=None, seed=10):
        self.seq2fitness = seq2fitness
        self.WT = WT
        self.AA_options = AA_options
        self.num_mut = num_mut
        self.mut_rate = mut_rate
        self.num_restarts = num_restarts
        self.max_steps = max_steps
        self.seq2fitness_many = seq2fitness_many
        random.seed(seed)

    def optimize(self):
        if self.seq2fitness_many is None or self.seq2fitness_many([self.WT] * 5) is None:
            def seq2fitness_many(seqs):
                return [self.seq2fitness(seq) for seq in seqs]

            self.seq2fitness_many = seq2fitness_many

        fit = self.seq2fitness(self.WT)
        best_best_muts = [self.WT[0] + '0' + self.WT[0], fit]
        self.fitness_trajectory = []

        for restart in range(
                self.num_restarts):  # At each restart, select a random set of num_mut mutations to begin from
            print('Beginning restart ' + str(restart))
            start_mut = generate_random_mut(self.WT, self.AA_options, self.num_mut)
            best_muts = [start_mut, self.seq2fitness(mut2seq(self.WT, start_mut))]
            best_fitnesses = []

            for step in range(self.max_steps):  # Set max number of steps in hill-climbing
                neighbors = []
                current_muts = best_muts[0]
                best_fitnesses.append(self.seq2fitness(mut2seq(self.WT, current_muts)))
                # Find all neighbors (where one mutations has been deleted and one added)
                for mut in current_muts.split(','):
                    for mut_pos in range(len(self.WT)):
                        if mut_pos not in [mut[1:-1] for mut in current_muts.split(',')]:
                            AA_options_pos = [AA for AA in self.AA_options[mut_pos] if AA is not self.WT[mut_pos]]
                            new_neighbor_muts = [self.WT[mut_pos] + str(mut_pos) + AA for AA in AA_options_pos]
                            new_neighbors = [current_muts.replace(mut, new_mut) for new_mut in new_neighbor_muts]
                            neighbors.extend(new_neighbors)
                # Predict the fitnesses for all neighbors
                pred_functions = list(self.seq2fitness_many([mut2seq(self.WT, muts) for muts in neighbors]))
                best_ind = pred_functions.index(max(pred_functions))
                print('Current fitness: ' + str(round(pred_functions[best_ind], 4)))
                if pred_functions[best_ind] == best_muts[1] or neighbors[best_ind] == best_muts[0]:
                    break
                best_muts = [neighbors[best_ind], pred_functions[best_ind]]
            self.fitness_trajectory.append(best_fitnesses)
            if best_muts[1] > best_best_muts[1]:
                best_best_muts = best_muts.copy()
        self.best_seq = [best_best_muts[0], best_best_muts[1]]
        return self.best_seq

    def plot_trajectory(self, savefig_name=None):
        for traj in self.fitness_trajectory:
            plt.plot(traj)
        plt.xlabel('Step')
        plt.ylabel('Fitness')
        plt.legend(['Restart ' + str(i) for i in range(self.num_restarts)])
        if savefig_name is None:
            plt.show()
            plt.close()
        else:
            plt.savefig(savefig_name)