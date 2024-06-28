import os

import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import cdist
from sklearn.cluster import DBSCAN, AgglomerativeClustering
import  numpy as np

from src.utils import av_gfp_WT


def practice_with_dendrograms():
    '''
    practice function for designing and making dendrograms with guasian distribution.
    :return:
    '''
    X,y= load_guassian_data()
    fig,ax= plt.subplots()

    colors= ['r','g','b']
    for i in np.arange(3):
        X_filter= X[i == y,:]
        ax.scatter(X_filter[:,0],X_filter[:,1],alpha=0.3,label=f'{i}',c=colors[i])

    ax.set_title('3 guassian distribution')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    fig.legend()
    fig.savefig(os.path.join('analysis','unit_test','3_guassian_dist.png'))


    dist_matrix =cdist(X,X,metric='euclidean')

    clustering_object= AgglomerativeClustering(n_clusters=3, metric='precomputed',
                                               linkage='average', compute_distances= True)

    predicted_y=clustering_object.fit_predict(dist_matrix)

    fig, ax = plt.subplots()
    colors = ['r', 'g', 'b']
    for i in np.arange(3):
        X_filter = X[i == predicted_y, :]
        ax.scatter(X_filter[:, 0], X_filter[:, 1], alpha=0.3, label=f'{i}', c=colors[i])
    ax.set_title('3 guassian distribution-- predicted labels')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    fig.legend()
    fig.savefig(os.path.join('analysis', 'unit_test', '3_guassian_dist_predicted.png'))

    fig,ax=plt.subplots()
    plot_dendrogram(clustering_object, ax=ax)
    ax.set_title('dendrogram')
    fig.savefig(os.path.join('analysis', 'unit_test', '3_guassian_dist_dendrogram.png'))
def load_guassian_data():
    '''
    generate fake guassian data to practice with dendrogram.
    :return: all the samples and corresponding labels
    '''
    # Set the mean vectors and covariance matrices for each 2D Gaussian distribution
    mean1 = np.array([0, 0])
    cov1 = np.array([[1, 0], [0, 1]])

    mean2 = np.array([10, 10])
    cov2 = np.array([[2, 0], [0, 2]])

    mean3 = np.array([-5, -5])
    cov3 = np.array([[0.5, 0], [0, 0.5]])

    # Set the number of samples to generate
    num_samples = 1000

    # Generate samples from the three 2D Gaussian distributions
    samples1 = np.random.multivariate_normal(mean1, cov1, num_samples)
    samples2 = np.random.multivariate_normal(mean2, cov2, num_samples)
    samples3 = np.random.multivariate_normal(mean3, cov3, num_samples)


    labels =np.array( [0]*len(samples1) + [1]*len(samples2) + [2]*len(samples3))
    # Concatenate the samples from all three distributions
    all_samples = np.concatenate((samples1, samples2, samples3))
    return all_samples,labels


def dendrogram_clustering_practice():
    '''
    fake data to practice clustering the dendrogram.
    :return:
    '''
    ytdist = np.array([662., 877., 255., 412., 996., 295., 468., 268.,
                       400., 754., 564., 138., 219., 869., 669.])
    Z = hierarchy.linkage(ytdist, 'single')
    plt.figure()
    dn = hierarchy.dendrogram(Z)
    plt.show()
def dbscan_clustering_precomputed(matrix_path,save_path,eps=10,min_samples=10):
    '''
    loading in precomputer clustering
    :param matrix_path: path to matrix if dbscan was precomputed
    :param save_path: file
    :param eps: eps value for sklearn DBscan
    :param min_samples: minimum number of samples in DBscan
    :return:
    '''
    dist_matrix = np.load(matrix_path)
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = dbscan.fit_predict(dist_matrix)


def my_own_test():
    '''
    verify that all the variants we are going to
     test experimentally have the correct mutations
    :return:
    '''
    df= pd.read_csv(os.path.join('results','final_chosen_sequences','final_sequences_test.csv'))

    for row in df.itertuples():
        if type(row.best_mutant)!=str:
            assert av_gfp_WT==row.best_sequence,'wt must be wt'
            print('verified wildtype variant')
        else:
            variant =row.best_mutant
            for mutant in variant.split(','):
               pos,mut= int(mutant[1:-1]),mutant[-1]
               assert row.best_sequence[pos]==mut
            print('verified variant : ',row.best_mutant)
if __name__ == '__main__':
    my_own_test()




