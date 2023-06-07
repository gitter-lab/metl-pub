
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
from sklearn.cluster import DBSCAN
import time
from Bio.Align import substitution_matrices
from Bio import Align
from sklearn.cluster import AgglomerativeClustering
from sklearn_extra.cluster import KMedoids

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import AffinityPropagation

# from general_stats import sort_variant
BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM62_min= np.min(BLOSUM62)
BLOSUM62_max= np.max(BLOSUM62)
BLOSUM62_AA_libary=  list(BLOSUM62.alphabet)
BLOSUMAA2Idx_MAPPING = {c: i for i, c in enumerate(BLOSUM62_AA_libary)}
blosum62 = substitution_matrices.load("BLOSUM62")

# create a pure numpy blosum62 matrix
# this is necessary for the vectorized function which
blosum62_arr = np.zeros((len(blosum62.alphabet), len(blosum62.alphabet)))
for i, aa1 in enumerate(blosum62.alphabet):
    for j, aa2 in enumerate(blosum62.alphabet):
        blosum62_arr[i, j] = blosum62[aa1, aa2]

# create a map from amino acid --> index in blosum matrix
aa_to_int_map = {aa: idx for idx, aa in enumerate(blosum62.alphabet)}


# a function to convert aa sequences to integer sequences
def aa_seq_to_int_array(seq):
    return np.array([aa_to_int_map[aa] for aa in seq], dtype=int)
# blosum similarity, vectorized with numpy
def vectorized_blosum_similarity(seq1, seq2, blosum_matrix):
    return blosum_matrix[np.array(seq1,dtype=int), np.array(seq2,dtype=int)].sum()

def onehot2index(onehot_encoding):
    idxs= []
    for row in onehot_encoding.reshape(-1,20):
        idxs.append(np.argmax(row))
    return np.array(idxs)

def blosum62_seq_similarity(seq1_onehot,seq2_onehot):
    seq1, seq2 = seq1_onehot.reshape(-1,len(BLOSUM62_AA_libary)), seq2_onehot.reshape(-1,len(BLOSUM62_AA_libary)).T
    score=np.sum(seq1 *( BLOSUM62 @ seq2).T)
    return score

def blosum62_seq_similarity_strings(seq1,seq2):
    tot =0
    for aa1,aa2 in zip(seq1,seq2):
        tot+= BLOSUM62[aa1,aa2]
    return tot

def blosum_onehot(sequence):
    '''
    different alphabet for the blosum62 matrix
    :param sequence:
    :return:
    '''
    X = np.zeros((len(sequence), len(BLOSUM62_AA_libary)))
    for i,s in enumerate(sequence):
        X[i][BLOSUMAA2Idx_MAPPING[s]] = 1
    return X



def unit_test_blosum62_seq_similarity(top_sequences_path,save_path):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: blosum_onehot(x).flatten())
    q= BLOSUM62_AA_libary
    p  = av_gfp_WT

    seq1=  df.iloc[0]['best_sequence']
    seq2 = df.iloc[1]['best_sequence']
    onehot1= blosum_onehot(seq1).flatten()
    onehot2 = blosum_onehot(seq2).flatten()

    print(blosum62_seq_similarity(onehot1,onehot2))
    print(blosum62_seq_similarity_strings(seq1,seq2))

    seq1_int = aa_seq_to_int_array(seq1)
    seq2_int = aa_seq_to_int_array(seq2)

    print('sam: ',vectorized_blosum_similarity(seq1_int, seq2_int, blosum62_arr))

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = 'local'
    aligner.gap_score = -np.inf

    A ,B,At,Bt,Ct= [],[],[],[],[]
    for i in tqdm(np.arange(len(df))):
        seq1 = df.iloc[0]['best_sequence']
        seq2 = df.iloc[i]['best_sequence']

        seq1_int = aa_seq_to_int_array(seq1)
        seq2_int = aa_seq_to_int_array(seq2)

        seq1_onehot=  df.iloc[0]['onehot']
        seq2_onehot = df.iloc[i]['onehot']
        start= time.time()
        aligner.align(seq1,seq2)
        middle1 =time.time()
        val1= blosum62_seq_similarity(seq1_onehot,seq2_onehot)
        middle2 =time.time()
        val2= vectorized_blosum_similarity(seq1_int,seq2_int,blosum62_arr)
        middle3 = time.time()
        val3= blosum62_seq_similarity_strings(seq1,seq2)
        end=time.time()

        assert val3==val2==val1,'all these values must be equal'
        At.append(middle2-middle1)
        Bt.append(middle3-middle2)
        Ct.append(end-middle3)



    fig,ax= plt.subplots()

    ax.scatter(A,B,s=0.1)
    ax.set_xlabel('alignment score')
    ax.set_title('seq similarity score vs alignment score-- GAP \n'
                 'comparison to first sequence')
    ax.set_ylabel('seq similarity')

    fig.savefig(os.path.join(save_path,'alignment_vs_seq_similarity_gap.png'))


def clustering_precomputed(matrix_path,top_sequences_path,clustering_object,
                                       save_path,n_clusters=5,nb_points=None,
                           down_sample=None,unique_down_sample=None):
    df = pd.read_csv(top_sequences_path)
    dist_matrix = np.load(matrix_path)
    if nb_points is not None:
        dist_matrix = dist_matrix[:nb_points, :nb_points]
        df=df[:nb_points]

    normalized_dist_matrix = 1 - (dist_matrix - BLOSUM62_min * 273) / (BLOSUM62_max * 273 - BLOSUM62_min * 273)

    # clusters = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
    labels = clustering_object.fit_predict(normalized_dist_matrix)
    assert n_clusters >= np.unique(labels).shape[0]
    fig,ax= plt.subplots()
    ax.hist(np.sum(normalized_dist_matrix,axis=0)/ normalized_dist_matrix.shape[0],bins=100,alpha=0.5)
    fig.savefig(os.path.join(save_path,'dist_matrix_full_distribution.png'))


    analyze_clusters(df=df, X=normalized_dist_matrix, labels=labels, n_clusters=n_clusters, save_path=save_path,
                     isX_dist_matrix=True)

    #down sample
    if down_sample is not None:
        cluster_idx=np.array(pd.Series(labels).value_counts().index[:down_sample])
        filter=np.isin(labels,cluster_idx )
        df_downsampled= df[filter]
        labels_downsampled= labels[filter]
        normalized_dist_matrix_down_sample =normalized_dist_matrix[filter][:,filter]
        analyze_clusters(df=df_downsampled, X=normalized_dist_matrix_down_sample,
                         labels=labels_downsampled, n_clusters=down_sample,
                         save_path=os.path.join(save_path,f'down_sample_{down_sample}'),
                         isX_dist_matrix=True,
                         cluster_indexes=cluster_idx)

    if unique_down_sample is not None:
        cluster_idx = np.array(pd.Series(labels).value_counts().index[:down_sample])












def dist_matrix_analysis(matrix_path,save_path):
    dist_matrix = np.load(matrix_path)
    fig, ax = plt.subplots()
    dist2everypoint = dist_matrix.sum(axis=1) / (dist_matrix.shape[0] - 1)
    ax.hist(dist2everypoint, bins=100, alpha=0.4)
    ax.set_title('full dataset average of distance \n'
                 'from each point to every other')
    ax.set_xlabel('average distance')
    fig.savefig(os.path.join(save_path, 'blosum62_encoding_avg_distance.png'))
    fig, ax = plt.subplots()
    diagonal_distances = np.diagonal(dist_matrix)
    ax.hist(diagonal_distances, bins=100, alpha=0.5)
    ax.set_title('full dataset diagonal distribution')
    ax.set_xlabel('diagonal distances')
    fig.savefig(os.path.join(save_path, 'blosum62_diagonal_distances.png'))


def compute_distances(top_sequences_path,save_path):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['blosum_integer'] = df['best_sequence'].apply(lambda x: aa_seq_to_int_array(x))
    X = np.vstack(df['blosum_integer'].to_numpy())
    dist_matrix = pairwise_distances(X,metric=lambda u,v: vectorized_blosum_similarity(u,v,blosum_matrix=blosum62_arr))
    np.save(os.path.join(save_path,'vectorized_blosum_similarity_distance_matrix.npy'),dist_matrix)




def dbscan_clustering_precomputed(matrix_path,save_path,eps=10,min_samples=10):
    dist_matrix = np.load(matrix_path)
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = dbscan.fit_predict(dist_matrix)



def dbscan_clustering_one(top_sequences_path,save_path,eps=90,min_samples=5):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())
    df['aa_idx'] = df['onehot'].apply(lambda x: onehot2index(x))
    X = np.vstack(df['aa_idx'].to_numpy())




def kmeans_clustering_onehot(top_sequences_path,save_path,n_clusters=5,random_state=0,
                             down_sample=None):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())

    X = np.vstack(df['onehot'].to_numpy())

    #1 - complete the kmeans clustering
    kmeans= KMeans(n_clusters=n_clusters,random_state=random_state).fit(X)
    labels = kmeans.labels_
    print('kmeans done')
    analyze_clusters(df=df,X=X,labels=labels,n_clusters=n_clusters,
                     isX_dist_matrix=False,
                     save_path=save_path)



    if down_sample is not None:
        cluster_idx=np.array(pd.Series(labels).value_counts().index[:down_sample])
        filter=np.isin(labels,cluster_idx)
        df_downsampled= df[filter]
        labels_downsampled= labels[filter]
        X_down_sample =X[filter,:]
        analyze_clusters(df=df_downsampled, X=X_down_sample,
                         labels=labels_downsampled, n_clusters=down_sample,
                         save_path=os.path.join(save_path,f'down_sample_{down_sample}'),
                         isX_dist_matrix=False,
                         cluster_indexes=cluster_idx)


def analyze_clusters(df,X,labels,n_clusters,save_path,isX_dist_matrix=True,distance_metric='hamming',
                     cluster_indexes=None):
    '''
    analyze a given set of clusters, points as onehot encodings
    :param df:
    :param X:
    :param labels:
    :param n_clusters:
    :param save_path:
    :return:
    '''


    training_variants = preprocess_train_variants_for_sim_anneal()
    arr=labels
    if cluster_indexes is None:
        cluster_indexes= np.arange(n_clusters)
    stats_df = pd.DataFrame(columns=['cluster', 'nb_points', 'inertia', 'min_variant', 'in_training_set',
                                     'train_variants'],
                            index=cluster_indexes)
    for nb in tqdm(cluster_indexes):

        filtered_df = df[arr == nb]

        # dist_matrix = cdist(filtered_X, filtered_X, metric='euclidean')
        if not isX_dist_matrix:
            filtered_X = X[arr == nb]
            dist_matrix = cdist(filtered_X, filtered_X, metric='euclidean')
        else:
            dist_matrix=X[arr==nb][:,arr==nb]

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

        normalized_sum_dist = sum_dist / (len(filtered_df))  # subtract one as to not include onself.
        stats_df.loc[nb] = [nb, len(sum_dist), np.sum(normalized_sum_dist), minimum_variant, count, tv]

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
    ax2.scatter(np.arange(len(cluster_indexes)),stats_df['in_training_set'],color='red', marker='o',label='in_training_set')
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

    #
    # runs2include = [run_params[0]]
    # for run_name in runs2include:
    #     fast_dbscan_clustering(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    #                                                           'analysis','top_sequences.csv'),
    #                           save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
    #                                                  'dbscan'),
    #                           )



    # compute_distance_matrix(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    #                                                           'analysis','top_sequences.csv'),
    #                         save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    #                                                           'analysis'))

    run_params = ['5']

    for run_name in tqdm(run_params):
        # compute_distances(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
        #                                                           'analysis','top_sequences.csv'),
        #                                   save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
        #                                                                        'analysis'))
        # dist_matrix_analysis(matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
        #                                                           'analysis','vectorized_blosum_similarity_distance_matrix.npy'), save_path=
        #                                     os.path.join('results', f'3d_{run_name}_mutant_10k_run','analysis'))
                     # )

        # eps = 0.4365
        # linkage='complete'
        # suffix= linkage
        # method = 'agglomerate'
        n_clusters = 50
        down_sample = 5
        save_dir = 'kmedioids'
        # clustering_object= AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage=linkage)
        # clustering_object= DBSCAN(eps=eps,metric='precomputed')
        clustering_object  = KMedoids(n_clusters=n_clusters,metric='precomputed',random_state=0)




        clustering_precomputed(clustering_object=clustering_object,
            matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                                  'analysis','vectorized_blosum_similarity_distance_matrix.npy'),
                                  save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
                                                         'analysis',save_dir,f'blosum62_clusters_{n_clusters}'),
                                           top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                              'analysis','top_sequences.csv'),
                               n_clusters=n_clusters,
                               down_sample=down_sample)




        # save_dir='kmeans'
        # kmeans_clustering_onehot(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
        #                                                       'analysis','top_sequences.csv'),
        #             save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
        #                                                  'analysis',save_dir,f'euclidean_clusters_{n_clusters}'),
        #             n_clusters=n_clusters,
        #             down_sample=down_sample,
        #             random_state=0)




    # unit_test_blosum62_seq_similarity(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    #                                                           'analysis','top_sequences.csv'),
    #                                   save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    #                                                                       'analysis'))
    # agglomerate_clustering_precomputed(matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    #                                                           'analysis','normalized_blosum62_dist_matrix.npy'),
    #                           save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
    #                                                  'analysis','agglomerate'))











    #
    # print(faster_blosum62_distance_normalized([4, 56, 8,4,76,9,7,40,5,4,79,7], [4, 56, 8,6,89,10,7,79,8,6,70,8]))
    # print(faster_blosum62_distance_normalized([4, 56, 8, 4, 76, 9, 7, 40, 5, 4, 79, 7],
    #                                           [4, 56, 8,4,76,9,7,40,5,4,79,7])) # s
    # print(fast_blosum62_distance_normalized([4, 56, 8,4,76,9,7,40,5,4,79,7], [4, 56, 8,6,89,10,7,79,8,6,70,8]))
        # , save_path, eps=0.5, min_samples=5)
        # kmeans_clustering_onehot(top_sequences_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
        #                                                          'analysis','top_sequences.csv'),
        #                          save_path=os.path.join('results',
        #                                                 f'3d_{run_name}_mutant_10k_run','analysis'))


