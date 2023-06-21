
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
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram
from general_stats import sort_variant,unique_sequences
import sys
sys.setrecursionlimit(10000)

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

def calculated_dist_matrix_onehot(top_sequences_path,save_path):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())
    X = np.vstack(df['onehot'].to_numpy())
    start= time.time()
    dist_matrix_onehot = cdist(X, X, metric='hamming')
    stop  = time.time()
    print(stop-start)
    np.save(os.path.join(save_path,'onehot_dist_matrix.npy'),dist_matrix_onehot)


def calculate_normalized_dist_matrix(dist_matrix):
    WT_BLOSUM = blosum62_seq_similarity_strings(
        'SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK',
        'SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK')
    calculated_max = (BLOSUM62_max - 4) * 10 + WT_BLOSUM
    calculated_min = 0
    normalized_dist_matrix = 1 - ((dist_matrix - calculated_min) / (calculated_max - calculated_min))
    return normalized_dist_matrix
def similarity_matrix_histograms(matrix_path,onehot_matrix_path,top_sequences_path,save_path):
    df = pd.read_csv(top_sequences_path)
    n_mutants = len(df.iloc[0]['best_mutant'].split(','))
    dist_matrix = np.load(matrix_path)

    normalized_dist_matrix= calculate_normalized_dist_matrix(dist_matrix)

    fig, ax = plt.subplots()
    ax.hist(normalized_dist_matrix.reshape(-1), bins=200, alpha=0.5)
    ax.set_title('distance matrix histogram')
    ax.set_xlabel('distance')
    ax.set_ylabel('count')
    fig.savefig(os.path.join(save_path, 'dist_matrix_full_distribution.png'))


    onehot_dist_matrix = np.load(onehot_matrix_path)*len(av_gfp_WT)*20
    fig, ax = plt.subplots()
    ax.hist(onehot_dist_matrix.reshape(-1), bins=20, alpha=0.5)
    ax.set_title('distance matrix histogram hamming distance full distribution')
    ax.set_xlabel('distance')
    ax.set_ylabel('count')
    fig.savefig(os.path.join(save_path, 'onehot_dist_matrix_full_distribution.png'))







def clustering_precomputed(matrix_path,
                           top_sequences_path,
                           onehot_matrix_path,
                           clustering_object,
                           save_path,
                           n_clusters=5,
                           down_sample=None,
                           percentage2include=0.15,
                           greedy_down_sample= None,
                           threshold=None,
                           n= 5000):

  


    df = pd.read_csv(top_sequences_path)[:n]
    df['best_mutant'] = df['best_mutant'].apply(lambda x:sort_variant(x))
    dist_matrix = np.load(matrix_path)[:n,:n]
    onehot_dist_matrix = np.load(onehot_matrix_path)*20*len(av_gfp_WT)
    onehot_dist_matrix=onehot_dist_matrix[:n,:n]

    print('length of dataframe: ',len(df))
    normalized_dist_matrix  = calculate_normalized_dist_matrix(dist_matrix)

    labels=clustering_object.fit_predict(normalized_dist_matrix)

    # labels= clustering_object.predict(normalized_dist_matrix)
    assert n_clusters >= np.unique(labels).shape[0]



    representative_sequences,mutations_df=analyze_clusters(df=df, global_dist_matrix=normalized_dist_matrix,
                     global_onehot_dist_matrix=onehot_dist_matrix,
                     labels=labels, n_clusters=n_clusters, save_path=save_path,
                     percentage2include=percentage2include)

    fig,ax=plt.subplots(figsize=(24,16))

    # plot_dendrogram(clustering_object, ax=ax)
    ax.set_title(f'n:{n} \n'
                 f' {save_path}',fontsize=30)
    fig.savefig(os.path.join(save_path,'dendrogram_full.png'))

    #down sample
    print('doing standard down sample')
    if down_sample is not None:
        cluster_idx=np.array(pd.Series(labels).value_counts().index[:down_sample])


        # filter=np.isin(labels,cluster_idx )
        # df_downsampled= df[filter]
        # labels_downsampled= labels[filter]
        # normalized_dist_matrix_down_sample =normalized_dist_matrix[filter][:,filter]
        analyze_clusters(df=df.copy(), global_dist_matrix=normalized_dist_matrix.copy(),
                         labels=labels.copy(), n_clusters=down_sample,
                         save_path=os.path.join(save_path,f'down_sample_{down_sample}'),
                         global_onehot_dist_matrix=onehot_dist_matrix.copy(),
                         cluster_indexes=cluster_idx)
    print('doing greedy down sample')
    if greedy_down_sample is not None:
        assert threshold is not None ,'threshold must be a positive integer '
        greedy_sequence_choice_algorithm(representative_variants=representative_sequences,
                                         global_dist_matrix=normalized_dist_matrix,
                                         df=df,
                                         save_path=os.path.join(save_path,f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}'),
                                         labels=labels,
                                         threshold=threshold,
                                         greedy_down_sample=greedy_down_sample)


    # if unique_down_sample is not None:
    #     cluster_idx = np.array(pd.Series(labels).value_counts().index[:down_sample])


def dendrogram_clustering_practice():
    ytdist = np.array([662., 877., 255., 412., 996., 295., 468., 268.,
                       400., 754., 564., 138., 219., 869., 669.])
    Z = hierarchy.linkage(ytdist, 'single')
    plt.figure()
    dn = hierarchy.dendrogram(Z)
    plt.show()








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




def kmeans_clustering_onehot(top_sequences_path,save_path,onehot_matrix_path,
                             n_clusters=5,random_state=0,
                             down_sample=None,percentage2include=0.15,
                            threshold =100, greedy_down_sample=None
                             ):

    df = pd.read_csv(top_sequences_path)
    df['best_mutant'] = df['best_mutant'].apply(lambda x: sort_variant(x))
    onehot_dist_matrix = np.load(onehot_matrix_path) * 20 * len(av_gfp_WT)
    onehot_dist_matrix = onehot_dist_matrix
    global_dist_matrix = np.sqrt(onehot_dist_matrix)

    fig, ax= plt.subplots()

    ax.hist(onehot_dist_matrix.reshape(-1),bins=20,alpha=0.5)
    ax.set_title('onehot full distribution')
    fig.savefig(os.path.join(save_path,'onehot_full_dist.png'))


    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())

    X = np.vstack(df['onehot'].to_numpy())

    #1 - complete the kmeans clustering
    kmeans= KMeans(n_clusters=n_clusters,random_state=random_state).fit(X)
    labels = kmeans.labels_
    print('kmeans done')

    representative_sequences, mutations_df= analyze_clusters(df=df,X=X,labels=labels,n_clusters=n_clusters,
                     save_path=save_path,
                     global_dist_matrix=global_dist_matrix,
                     global_onehot_dist_matrix=onehot_dist_matrix,
                     percentage2include=percentage2include)


    print('doing greedy down sample')
    if greedy_down_sample is not None:
        assert threshold is not None, 'threshold must be a positive integer '
        greedy_sequence_choice_algorithm(representative_variants=representative_sequences,
                                         global_dist_matrix=global_dist_matrix,
                                         df=df,
                                         save_path=os.path.join(save_path,
                                                                f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}'),
                                         labels=labels,
                                         threshold=threshold,
                                         greedy_down_sample=greedy_down_sample)

    # if down_sample is not None:
    #     cluster_idx=np.array(pd.Series(labels).value_counts().index[:down_sample])
    #     filter=np.isin(labels,cluster_idx)
    #     df_downsampled= df[filter]
    #     labels_downsampled= labels[filter]
    #     X_down_sample =X[filter,:]
    #     analyze_clusters(df=df_downsampled, X=X_down_sample,
    #                      labels=labels_downsampled, n_clusters=down_sample,
    #                      save_path=os.path.join(save_path,f'down_sample_{down_sample}'),
    #                      isX_dist_matrix=False,
    #                      cluster_indexes=cluster_idx)

def analyze_clusters(df,labels,n_clusters,save_path,
                     X=None,
                     global_dist_matrix=None,
                     global_onehot_dist_matrix=None,
                     cluster_indexes=None,
                     percentage2include=0.15):
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


    mutations_df= pd.DataFrame(index=cluster_indexes)
    representative_sequences = pd.DataFrame(index=cluster_indexes)
    for nb in tqdm(cluster_indexes):

        filtered_df = df[arr == nb].copy()
        filtered_df.to_csv(os.path.join(save_path,f'cluster_{nb}.csv'))

        unique_sequences(filtered_df['best_mutant'], other_unique=2, save_dir=save_path, suffix=f"{nb}_variants")

        mutations =','.join(list(filtered_df['best_mutant'])).split(',')
        unique_sequences(filtered_df['best_mutant'],
                            other_unique=2,
                            save_dir=save_path,
                            suffix=f"{nb}_variants")
        mutation_counts= unique_sequences(pd.Series(mutations),
                            other_unique=2,
                            save_dir=save_path,
                            suffix=f"{nb}_mutations")

        mutation_counts=mutation_counts[mutation_counts>len(filtered_df)*percentage2include*5]


        for mutation_name,mutation_count in zip(mutation_counts.index,mutation_counts):
            mutations_df.loc[nb,mutation_name]=int(mutation_count)





        # dist_matrix = cdist(filtered_X, filtered_X, metric='euclidean')
        if X is not None:
            filtered_X = X[arr == nb]
            dist_matrix = cdist(filtered_X, filtered_X, metric='euclidean')
        else:
            dist_matrix=global_dist_matrix[arr==nb][:,arr==nb]


        if global_onehot_dist_matrix is not None:
            onehot_dist_matrix =  global_onehot_dist_matrix[arr==nb][:,arr==nb]
            fig, ax = plt.subplots()
            ax.hist(onehot_dist_matrix.reshape(-1), bins=20, alpha=0.5)
            ax.set_title('distance matrix histogram \n '
                         'onehot distance full distribution')
            ax.set_xlabel('distance')
            ax.set_ylabel('count')
            fig.savefig(os.path.join(save_path, f'onehot_cluster_{nb}.png'))



        # sum up the pairwise distances
        sum_dist = np.sum(dist_matrix, axis=0)

        # agrument of the smallest point
        argmin = np.argmin(np.abs(sum_dist))

        # get the smallest variant
        minimum_variant = filtered_df.iloc[argmin]['best_mutant']

        representative_sequences.loc[nb,'best_mutant'] = minimum_variant
        representative_sequences.loc[nb, 'best_fitness'] = filtered_df.iloc[argmin]['best_fitness']


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

        ax.set_xlabel("pairwise distance every point (N x N)")
        ax.set_ylabel("frequency")

        # Set number of bins for histogram
        num_bins = 200

        # Plot histogram of data with specified number of bins

        n, bins, patches = ax.hist(dist_matrix.reshape(-1),
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
        ax.legend([f"Mininum (Average) Value = {min:.2f}"], loc="upper right")

        # Save plot as image file
        plt.tight_layout()
        fig.savefig(os.path.join(save_path, f'cluster_{nb}.png'), dpi=300)

    # nb_mutations=  len(df.iloc[0]['best_mutant'].split(','))
    # mutations_df = mutations_df/5
    mutations_df.to_csv(os.path.join(save_path,f'mutation_scaffold_percent_{percentage2include}.csv'))
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

    return representative_sequences,mutations_df
def mutual_information():
    pass
def silolette_factor():
    pass

def load_guassian_data():
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



def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    start= time.time()
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)
    print(time.time()-start)

def practice_with_dendrograms():
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


def greedy_sequence_choice_algorithm(representative_variants,
                                    global_dist_matrix,
                                     df,
                                     greedy_down_sample,
                                     save_path,
                                     labels,
                                     threshold):

    print('\n ===== starting greedy approach =====')
    filter_boolean= []
    # for column in mutations_df.columns:
    #     # remove the other column
    #     if 'other' in column:
    #         mutations_df= mutations_df.drop(column,axis=1)


    # only include clusters which are above a certain threshold

    # =representative_variants[]
    # for row in mutations_df.itertuples():
    #     filter_boolean.append(np.any(np.array(row) >=  threshold))


    #  update 6.20.23 include anything that has a cluster count above a threshold,
    #       not what is defined by the scaffold matrix
    representative_sequences= representative_variants[pd.Series(labels).value_counts() >=threshold]
    idx_of_rep_variants =[]
    for row in representative_sequences.itertuples():
        idx_of_rep_variants.append(np.argmax(df['best_mutant']==row.best_mutant))

    representative_sequences['idx']  =idx_of_rep_variants


    local_dist_matrix= pd.DataFrame(index=representative_sequences.index,
                                    columns=representative_sequences.index)

    for i,idxi in enumerate(representative_sequences.index):
        for j,idxj in enumerate(representative_sequences.index):
            local_dist_matrix.loc[idxi,idxj]=global_dist_matrix[idx_of_rep_variants[i],idx_of_rep_variants[j]]

    # filter down the distance matrix
    # local_dist_matrix_test = global_dist_matrix[np.array(idx_of_rep_variants)][:, np.array(idx_of_rep_variants)]
    print('stop')


    # find starting sequence

    chosen = []
    remaining = list(representative_sequences.index)
    for label in pd.Series(labels).value_counts().index:
        if int(label) in representative_sequences.index:
            chosen.append(int(label))
            remaining.remove(int(label))
            break


    print('chosen: ',chosen)
    print('remaining',remaining)

            # remove from remaining



    ## unit_test !!!!
    # local_dist_matrix  =pd.DataFrame(index=[4,6,8,10],columns=[4,6,8,10],data=[[0,9,7,4],[9,0,3,5],[7,3,0,2],[4,5,2,0]])
    # chosen = [4]
    # remaining=  [6,8,10]
    # greedy_down_sample = 3

    while len(chosen)<greedy_down_sample:
        avg_distance = []
        for r_point in remaining:
            tot_distance_point = 0
            for c_point in chosen:
                tot_distance_point += local_dist_matrix.loc[c_point,r_point]

            avg_distance_point = tot_distance_point / len(chosen)
            avg_distance.append(avg_distance_point)
            print('\n avg_dist : ',avg_distance)

        idx_max= np.argmax(avg_distance)

        variant2choose= remaining[idx_max]
        chosen.append(variant2choose)
        remaining.remove(variant2choose)

        print('\nchosen : ',chosen)
        print('remaining : ',remaining)


    chosen_sequences =representative_sequences.loc[chosen].reset_index(names='cluster')

    variants_list = []
    for variant in chosen_sequences['best_mutant']:
        variants_list.append(pd.Series([variant]))
    comparison_labels = list(chosen_sequences['cluster'])



    aa_distribution_comparison(variants_list,comparison_labels,save_path)
    res_distribution_comparison(variants_list,comparison_labels,save_path)


    mutations = ",".join(chosen_sequences['best_mutant']).split(',')
    residues = [int(mutation[1:-1]) for mutation in mutations]
    nb_of_unique_mutations = len(np.unique(mutations))
    nb_of_unique_residues = len(np.unique(residues))

    chosen_sequences['unique_mutations']= nb_of_unique_mutations
    chosen_sequences['unique_residues']= nb_of_unique_residues


    chosen_sequences.to_csv(os.path.join(save_path,'chosen_sequences.csv'))
    # initilize remaining representative clusters to be chosen...
        # while len(list_chosen) < nb_points2choose:
        # avg_distance = [ ]
        # for each point in remaining:
        # tot_distance_point = 0
        # for each point in list chosen :
        #
    # avg_distance_point = tot_distance_point / len(list_chosen)
    # avg_distance.append(avg_distance_point)




    # filter down the clusters using the cluster matrix ,  then pick a representative
    # sequence from each of the clusters by filtering down the global distance matrix
    # and then get filter down the global distance matrix to just the cluster representative
    # sequences and use that distance matrix in the algorithm

    # assert that the number of sequences that meet the criteria >= nb of sequences your hoping to get

    # choose the sequence with the highest number points in its cluster
    # for each representative point that is remaining:







    # save all the best mutants chosen

    # do amino acid analysis on them

    # do residue analysis on them

    # plot the distribution of score and cluster number for each selected one


    pass





def run_cluster_precomputed():
    run_params = ['rejection_5','only_train_5']
    Precentage2include = [0.05,0.15]
    n_clusters = [10,20,50]
    linkage ='complete'
    suffix = linkage
    method = 'agglomerate'
    n_reduce = 10000
    reduce_large_set = ''
    if n_reduce < 10000:
        reduce_large_set = f'_{n_reduce}'


    down_sample = None
    greedy_down_sample = 5
    save_dir = method + '_' + suffix
    distance = 'blosum62_clusters'
    threshold = 100




    for n_cluster in n_clusters:
        for run_name,percentage2include in zip(run_params,Precentage2include):
                #     # compute_distances(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                           'analysis','top_sequences.csv'),
                #     #                                   save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                                        'analysis'))
                #     # dist_matrix_analysis(matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                           'analysis','vectorized_blosum_similarity_distance_matrix.npy'), save_path=
                #     #                                     os.path.join('results', f'3d_{run_name}_mutant_10k_run','analysis'))
                #                  # )
                #
                # eps = 0.4365


                #     ## 1)
                #     # calculated_dist_matrix_onehot(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                       'analysis','top_sequences.csv'),
                #     #                                save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                  'analysis'))
                #
                #     # # 2)
                #     # similarity_matrix_histograms(matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                           'analysis','vectorized_blosum_similarity_distance_matrix.npy'),
                #     #                           save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
                #     #                                                  'analysis'),
                #     #                                    top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                #     #                                                       'analysis','top_sequences.csv'),
                #     #                              onehot_matrix_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
                #     #                                                  'analysis','onehot_dist_matrix.npy'))
                #
                #
                #
                clustering_object = AgglomerativeClustering(n_clusters=n_cluster,
                                                            metric='precomputed',
                                                            linkage=linkage,
                                                            compute_distances=True)
                #     # # clustering_object= DBSCAN(eps=eps,metric='precomputed')
                #     # clustering_object  = KMedoids(n_clusters=n_clusters,metric='precomputed',random_state=0)

                #     ## make a directory if necessary
                clustering_directory = os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                    f'analysis{reduce_large_set}', save_dir)
                if not os.path.exists(clustering_directory):
                    os.mkdir(clustering_directory)

                distance_directory = os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                  f'analysis{reduce_large_set}', save_dir, f'{distance}_{n_cluster}')
                if not os.path.exists(distance_directory):
                    os.mkdir(distance_directory)

                if not os.path.exists(os.path.join(distance_directory, 'down_sample_5')):
                    os.mkdir(os.path.join(distance_directory, 'down_sample_5'))

                if not os.path.exists(os.path.join(distance_directory, f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}')):
                    os.mkdir(os.path.join(distance_directory, f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}'))




                if 'onehot' in distance:
                    matrix2use = 'onehot_dist_matrix.npy'
                else:
                    matrix2use = 'vectorized_blosum_similarity_distance_matrix.npy'

                print('matrix: ', matrix2use)

                clustering_precomputed(
                    clustering_object=clustering_object,
                    matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                             'analysis', matrix2use),
                    save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                           f'analysis{reduce_large_set}', save_dir, f'{distance}_{n_cluster}'),
                    top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                    'analysis', 'top_sequences.csv'),
                    onehot_matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                    'analysis', 'onehot_dist_matrix.npy'),
                    n_clusters=n_cluster,
                    down_sample=down_sample,
                    percentage2include=percentage2include,
                    greedy_down_sample=greedy_down_sample,
                    threshold=threshold,
                    n  =n_reduce
                )

def run_kmeans_onehot():
    run_name = 'only_train_5'
    n_clusters = 5
    save_dir = 'kmeans'
    distance = 'euclidean'
    down_sample=None
    threshold= 100
    greedy_down_sample=None

    clustering_directory = os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                        'analysis', save_dir)
    if not os.path.exists(clustering_directory):
        os.mkdir(clustering_directory)

    distance_directory = os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                      'analysis', save_dir, f'{distance}_{n_clusters}')
    if not os.path.exists(distance_directory):
        os.mkdir(distance_directory)



    if not os.path.exists(
            os.path.join(distance_directory, f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}')):
        os.mkdir(os.path.join(distance_directory, f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}'))

    if not os.path.exists(os.path.join(distance_directory, 'down_sample_5')):
        os.mkdir(os.path.join(distance_directory, 'down_sample_5'))


    matrix2use = 'onehot_dist_matrix.npy'

    kmeans_clustering_onehot(
        save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                               'analysis', save_dir, f'{distance}_{n_clusters}'),
        top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                        'analysis', 'top_sequences.csv'),
        onehot_matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                        'analysis', 'onehot_dist_matrix.npy'),
        n_clusters=n_clusters,
        down_sample=down_sample,
        random_state=0,
        greedy_down_sample=greedy_down_sample,
        threshold=threshold,
    percentage2include=0.05)





if __name__ == '__main__':
    # run_kmeans_onehot()
    run_cluster_precomputed()
    # practice_with_dendrograms()
#   # print('stop')
    #
    #
    #
    #
    #
    # # dendrogram_clustering_practice()
    # # runs2include = [run_params[0]]
    # # for run_name in runs2include:
    # #     fast_dbscan_clustering(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    # #                                                           'analysis','top_sequences.csv'),
    # #                           save_path=os.path.join('results',f'3d_{run_name}_mutant_10k_run',
    # #                                                  'dbscan'),
    # #                           )
    #
    #
    #
    # # compute_distance_matrix(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    # #                                                           'analysis','top_sequences.csv'),
    # #                         save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
    # #                                                           'analysis'))
    #
    #
    #









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


