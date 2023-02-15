
### this script will run k-means from scikit-learn,
### also maybe where we can visualize where the sequences are
## if we run PCA! using embeddings for METL-L

import matplotlib.pyplot as plt
import numpy as np


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
def kmeans_clustering():
    pass

def dimensionality_reduction():
    '''
    plot clusters with dimensionality reduction
    :return:
    '''
    pass
def plot_trajetory_dimensionality_reduction():
    pass

def hamming_distance_comparison_with_datasets():
    '''
    see how much the hamming distance compares with the generated mutants and the sequences in the
    dataset.
    :return:
    '''
    pass


if __name__ == '__main__':
    pass

