
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from unit_tests.preprocessing_gfp_run import preprocess_train_variants_for_sim_anneal
from tqdm import tqdm
from comparison_stats import aa_distribution_comparison,res_distribution_comparison
from scipy.spatial.distance import cdist
import os
from sklearn.cluster import KMeans
from src.utils import onehot,variant2sequence,av_gfp_WT,onehot2sequence
from Bio.Align import substitution_matrices
from sklearn.cluster import DBSCAN
import time
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import AffinityPropagation
from general_stats import sort_variant
BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM62_min= np.min(BLOSUM62)
BLOSUM62_max= np.max(BLOSUM62)
BLOSUM62_AA_libary=  BLOSUM62.alphabet

# Define your custom distance metric

baseline_score =np.sum([BLOSUM62[x,x] for x in av_gfp_WT])


def unit_tests_distance_tests(top_sequences_path):
    df = pd.read_csv(top_sequences_path)
    df['blosum_encoding'] = df['best_mutant'].apply(lambda x: variant2fastblosum(x))
    print('baseline score', baseline_score)
    idx1,idx2= 0,1
    variant1 = df.iloc[idx1]['best_mutant']
    variant2=  df.iloc[idx2]['best_mutant']



    total = baseline_score
    for mutant in variant1.split(','):
        total+= BLOSUM62[mutant[0],mutant[-1]]  - BLOSUM62[mutant[0],mutant[0]]
    for mutant in variant2.split(','):
        total+= BLOSUM62[mutant[0],mutant[-1]] - BLOSUM62[mutant[0],mutant[0]]

    print('total score no overlap ',total)

    print('fast blosum score ',fast_blosum62_distance(df.iloc[idx1]['blosum_encoding'],df.iloc[idx2]['blosum_encoding']))




    print('=======positional =========')
    print('fast blosum score ',
          fast_blosum62_distance([4,56,8], [4,56,8]))

    print('positional check:',
          BLOSUM62[8,8] - BLOSUM62[4,4] + baseline_score)



    # check looks like everything checks out
def unit_test_normalized_blosum62(top_sequences_path):
    # df = pd.read_csv(top_sequences_path)
    # df['blosum_encoding'] = df['best_mutant'].apply(lambda x: variant2fastblosum(x))
    # print('baseline score', baseline_score)

    print('positional ')
    print('fast blosum score normalized',
          fast_blosum62_distance_normalized([4, 56, 8], [4, 56, 8]))

    print('positional check:',
          1- ((BLOSUM62[8, 8] - BLOSUM62_min)/ (BLOSUM62_max - BLOSUM62_min)))

    print('mutational')
    print('fast blosum score normalized ',
          fast_blosum62_distance_normalized([8, 57, 8], [8, 56, 8]))

    print('mutational check:',
          (1 - ((BLOSUM62[8, 8] - BLOSUM62_min) / (BLOSUM62_max - BLOSUM62_min)))*2)


def fast_blosum62_distance_normalized(seq1,seq2):
    '''
    this is the normalized
    :param seq1: blosum62 encoding which is the [wt_idx ,  position, mutation_idx, ... ] where index is specified by
    BLOSUM62_AA_libary
    :param seq2:  blosum62 encoding which is the [wt_idx ,  position, mutation_idx , ...  ] where index is specified by
    BLOSUM62_AA_libary
    :return: normalized blosum62 distance metric where lower is more similarity (max is 2*number of mutations)
    '''
    seq1_wt, seq1_pos, seq1_mut = seq1[::3], seq1[1::3], seq1[2::3]
    seq2_wt, seq2_pos, seq2_mut = seq2[::3], seq2[1::3], seq2[2::3]
    tot_score= 0

    for wt, pos, mut in zip(seq1_wt, seq1_pos, seq1_mut):
        indice_of_intersection = np.where(np.in1d(seq2_pos, [pos]))[0]
        if len(indice_of_intersection) > 0:
            # means we are at the same position [ mutation -> mutation]
            idx = indice_of_intersection[0]
            mut2 =  seq2_mut[idx]
            # 1 -  [0,1] where 1 is the most similar so we want that to be zero.
            tot_score += 1- (BLOSUM62[mut, mut2] - BLOSUM62_min) / (BLOSUM62_max - BLOSUM62_min)

            # remove this mutation so it is no longer included
            seq2_mut = np.delete(seq2_mut,idx)
            seq2_wt= np.delete(seq2_wt, idx)
            seq2_pos =np.delete(seq2_pos, idx)


            # seq2_mut.pop(idx)
            # seq2_wt.pop(idx)
            # seq2_pos.pop(idx)
        else:
            # means we are at a different position [ wildtype -> mutation]
            tot_score += 1-((BLOSUM62[mut, wt] - BLOSUM62_min) / (BLOSUM62_max - BLOSUM62_min))

        # no overlaping positions [wildtype -> mutation]
    for wt2, pos2, mut2 in zip(seq2_wt, seq2_pos, seq2_mut):
        tot_score += 1-((BLOSUM62[mut2, wt2] - BLOSUM62_min) / (BLOSUM62_max - BLOSUM62_min))

    return tot_score

def fast_blosum62_distance(seq1,seq2):
    '''
    fast blosum62 scoring metric
    representation given as [index aa wt, position, index aa mutant,
                                index aa wt,position, index aa mutant , ... ]
    :param seq1: representation as shown above
    :param seq2: representation as shown above
    :return: blosum62 substitution score for an entire sequence (specificially avgfp sequence)
    '''
    tot_score = baseline_score  # for example
    seq1_wt, seq1_pos, seq1_mut =seq1[::3] ,seq1[1::3], seq1[2::3]
    seq2_wt, seq2_pos, seq2_mut =seq2[::3] , seq2[1::3], seq2[2::3]


    for wt,pos,mut in zip(seq1_wt,seq1_pos,seq1_mut):
        indice_of_intersection= np.where(np.in1d(seq2_pos,[pos]))[0]
        if len(indice_of_intersection) > 0:
            # means we are at the same position [ mutation -> mutation]
            idx =indice_of_intersection[0]
            wt2,pos2, mut2= seq2_wt[idx],seq2_pos[idx], seq2_mut[idx]
            assert wt2==wt ,'these should be equal'
            assert pos ==pos2,'these should be equal'
            tot_score += BLOSUM62[mut,mut2] - BLOSUM62[wt,wt]

            # remove this mutation so it is no longer included
            seq2_mut.pop(idx)
            seq2_wt.pop(idx)
            seq2_pos.pop(idx)
        else:
            # means we are at a different position [ wildtype -> mutation]
            tot_score+= BLOSUM62[mut,wt] - BLOSUM62[wt,wt]

    # no overlaping positions [wildtype -> mutation]
    for wt2,pos2,mut2 in zip(seq2_wt,seq2_pos,seq2_mut):
        tot_score += BLOSUM62[mut2, wt2] - BLOSUM62[wt2, wt2]


    return tot_score

def blosum62_distance(seq1_onehot,seq2_onehot):
    tot_score = 0
    for aa1,aa2 in zip(seq1_onehot,seq2_onehot):
        blosum_score = BLOSUM62[aa1,aa2]
        normalized_score = 1 - (blosum_score - BLOSUM62_min ) / ( BLOSUM62_max- BLOSUM62_min)
        tot_score+= normalized_score
    return tot_score



def onehot2index(onehot_encoding):
    idxs= []
    for row in onehot_encoding.reshape(-1,20):
        idxs.append(np.argmax(row))
    return np.array(idxs)


def variant2fastblosum(variant):
    encoding =[ ]

    for mutant in variant.split(','):
        encoding.append(BLOSUM62_AA_libary.find(mutant[0]))
        encoding.append(int(mutant[1:-1]))
        encoding.append(BLOSUM62_AA_libary.find(mutant[-1]))
    return encoding



def compute_distance_matrix(top_sequences_path,save_name):
    length = 200
    df = pd.read_csv(top_sequences_path)[:length]
    df['blosum_encoding']  = df['best_mutant'].apply(lambda x: variant2fastblosum(x))
    X = np.vstack(df['blosum_encoding'].to_numpy())
    start= time.time()
    dist_matrix = cdist(X, X, metric=lambda u, v: fast_blosum62_distance_normalized(u, v))
    end =time.time()
    print('time : ',end-start)
    np.save(os.path.join('analysis','distance_matrix_results',f'{save_name}_normalized_blosum62_dist_matrix.npy'),
            dist_matrix)
    # clustering = DBSCAN(eps=eps, min_samples=min_samples, metric=lambda u, v: fast_blosum62_distance_normalized(u, v)).fit(X)
    # labels = clustering.labels_
    # print(labels)

def dbscan_clustering_one(top_sequences_path,save_path,eps=90,min_samples=5):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())
    df['aa_idx'] = df['onehot'].apply(lambda x: onehot2index(x))
    X = np.vstack(df['aa_idx'].to_numpy())
    dist_matrix = cdist(X, X, metric=lambda u,v : blosum62_distance(u,v))

    # sum up the pairwise distances
    sum_dist = np.sum(dist_matrix, axis=0)

    fig,ax = plt.subplots(1,1)

    ax.hist(sum_dist,bins=50,alpha=0.5)

    fig.savefig(os.path.join(save_path,'entire_distribution_blosumscore.png'))


    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric=lambda u,v : blosum62_distance(u,v)).fit(X)
    labels = clustering.labels_


    print('stop')


def kmeans_clustering_onehot(top_sequences_path,save_path,n_clusters=5,random_state=0):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())

    X = np.vstack(df['onehot'].to_numpy())

    #1 - complete the kmeans clustering
    kmeans= KMeans(n_clusters=n_clusters,random_state=random_state).fit(X)
    labels = kmeans.labels_

    analyze_clusters(df=df,X=X,labels=labels,n_clusters=n_clusters,
                     save_path=save_path)


def analyze_clusters(df,X,labels,n_clusters,save_path,distance_metric='euclidean'):
    '''
    analyze a given set of clusters, points as onehot encodings
    :param df:
    :param X:
    :param labels:
    :param n_clusters:
    :param save_path:
    :return:
    '''
    stats_df = pd.DataFrame(columns=['cluster', 'nb_points', 'inertia', 'min_variant', 'in_training_set',
                                     'train_variants'],
                            index=np.arange(n_clusters))

    training_variants = preprocess_train_variants_for_sim_anneal()
    arr=labels
    for nb in tqdm(np.arange(n_clusters)):

        filtered_df = df[arr == nb]
        filtered_X = X[arr == nb]
        # dist_matrix = cdist(filtered_X, filtered_X, metric='euclidean')
        dist_matrix = cdist(filtered_X, filtered_X, metric=distance_metric)

        # sum up the pairwise distances
        sum_dist = np.sum(dist_matrix, axis=0)

        # agrument of the smallest point
        argmin = np.argmin(np.abs(sum_dist))

        # get the smallest variant
        minimum_variant = filtered_df.iloc[argmin]['best_mutant']

        # find the similarity of this point to the training set.
        count = 0
        tv = []
        for v in minimum_variant.split(','):
            if v in training_variants:
                tv.append(v)
                count += 1

        normalized_sum_dist = sum_dist / (len(filtered_df) - 1)  # subtract one as to not include onself.
        stats_df.iloc[nb] = [nb, len(sum_dist), np.sum(normalized_sum_dist), minimum_variant, count, tv]

        # Create figure and axes objects
        fig, ax = plt.subplots()

        # Set plot title and axis labels
        ax.set_title(f"cluster: {nb}, nb_points:{len(sum_dist)} , \n"
                     f"inertia (total sum) : {np.sum(sum_dist)}\n,"
                     f"mininum point: {minimum_variant},\n"
                     f"in_training_set : {count},\n"
                     f"train_variants : {tv}")

        ax.set_xlabel("pairwise distance")
        ax.set_ylabel("frequency")

        # Set number of bins for histogram
        num_bins = 200

        # Plot histogram of data with specified number of bins

        n, bins, patches = ax.hist(normalized_sum_dist,
                                   bins=num_bins, color="#6495ED", alpha=0.5)

        # Customize appearance of histogram
        for i in range(len(patches)):
            # Set edge color and line width of each bar in histogram
            patches[i].set_edgecolor("#000080")
            patches[i].set_linewidth(0.5)

        # Add vertical line for mean value
        min = normalized_sum_dist[argmin]
        ax.axvline(min, color="#FFA500", linestyle="--", linewidth=1.5)

        # Add legend for mean value
        ax.legend([f"Mininum (Normalized) Value = {min:.2f}"], loc="upper right")

        # Save plot as image file
        plt.tight_layout()
        fig.savefig(os.path.join(save_path, f'cluster_{nb}.png'), dpi=300)

    mutations = ",".join(stats_df['min_variant']).split(',')
    residues = [int(mutation[1:-1]) for mutation in mutations]
    nb_of_unique_mutations = len(np.unique(mutations))
    nb_of_unique_residues = len(np.unique(residues))
    fig, ax = plt.subplots()
    ax = stats_df[['nb_points']].plot.bar(ax=ax, color=['blue'], alpha=0.5)
    ax.set_title('Cluster Statistics\n'
                 f'Unique Mutations: {nb_of_unique_mutations} / {len(mutations)}\n'
                 f'Unique Residues: {nb_of_unique_residues} / {len(mutations)}')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Number in Cluster')

    # Create twin axis for "in_training_set" data
    ax2 = ax.twinx()
    stats_df['in_training_set'].plot.line(ax=ax2, color='red', marker='o')
    ax2.set_ylabel('In Training Set')
    # Set legend for both axes
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2)
    plt.tight_layout()
    # fig,ax  =plt.subplots()
    # ax=stats_df[['nb_points','in_training_set','inertia']].plot.bar(ax=ax)
    # ax.set_title('Cluster Statistics')
    # ax.set_xlabel('Cluster')
    # ax.legend()
    plt.savefig(os.path.join(save_path, f'cluster_stats.png'))

    # lets look at where the residues are in those chosen sequences

    # do some preprocessing on the min variant
    variants_list = []
    for variant in stats_df['min_variant']:
        variants_list.append(pd.Series([variant]))
    labels = list(stats_df['cluster'])

    res_distribution_comparison(variants_list, labels, save_dir=save_path)
    aa_distribution_comparison(variants_list, labels, save_dir=save_path)
def mutual_information():
    pass
def silolette_factor():
    pass

if __name__ == '__main__':

    run_params = ['5','rejection_5','only_train_5', '10','rejection_10', 'only_train_10']
    runs2include = [run_params[0]]
    for run_name in runs2include:
        fast_dbscan_clustering(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                              'analysis','top_sequences.csv'),
                              save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
                                                     'dbscan'),
                              )



        # , save_path, eps=0.5, min_samples=5)
        # kmeans_clustering_onehot(top_sequences_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
        #                                                          'analysis','top_sequences.csv'),
        #                          save_path=os.path.join('results',
        #                                                 f'3d_{run_name}_mutant_10k_run','analysis'))


