
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tqdm import tqdm
import time,os
import sys
# increase  recursion limit for increasing the size of the dendrogram
sys.setrecursionlimit(10000)


# clustering imports
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import cdist


from src.utils import onehot,variant2sequence,av_gfp_WT
from preprocessing_gfp_run import preprocess_train_variants_for_sim_anneal
from comparison_stats import aa_distribution_comparison,res_distribution_comparison
from general_stats import sort_variant,unique_sequences

from distances import compute_distances,calculated_dist_matrix_onehot,calculate_normalized_dist_matrix


# plot dendrogram
from scipy.cluster.hierarchy import dendrogram


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
                           n= 5000,
                           include_dendrogram = False ):

    '''
    the code to run
    :param matrix_path:
    :param top_sequences_path:
    :param onehot_matrix_path:
    :param clustering_object:
    :param save_path:
    :param n_clusters:
    :param down_sample:
    :param percentage2include:
    :param greedy_down_sample: int, number of sequences to choose , must be less than or equal to
                                    the number of unique integers in labels array
    :param threshold: number of sequences that must be in a clsuter to be considered
    :param n: number of top sequences to include
    :param include_dendrogram:boolean to make dendrogram
    :return:
    '''

    df = pd.read_csv(top_sequences_path)[:n]
    df['best_mutant'] = df['best_mutant'].apply(lambda x:sort_variant(x))
    dist_matrix = np.load(matrix_path)[:n,:n]
    onehot_dist_matrix = np.load(onehot_matrix_path)*20*len(av_gfp_WT)
    onehot_dist_matrix=onehot_dist_matrix[:n,:n]

    print('length of dataframe: ',len(df))
    normalized_dist_matrix  = calculate_normalized_dist_matrix(dist_matrix)

    labels=clustering_object.fit_predict(normalized_dist_matrix)

    assert n_clusters >= np.unique(labels).shape[0]



    representative_sequences,mutations_df=analyze_clusters(df=df, global_dist_matrix=normalized_dist_matrix,
                     global_onehot_dist_matrix=onehot_dist_matrix,
                     labels=labels, n_clusters=n_clusters, save_path=save_path,
                     percentage2include=percentage2include)

    if include_dendrogram:
        fig,ax=plt.subplots(figsize=(24,16))


        plot_dendrogram(clustering_object, ax=ax)
        ax.set_title(f'n:{n} \n'
                     f' {save_path}',fontsize=30)
        fig.savefig(os.path.join(save_path,'dendrogram_full.png'))

    #down sample
    print('doing standard down sample')
    if down_sample is not None:
        # down sample based on 5 most populous clusters
        cluster_idx=np.array(pd.Series(labels).value_counts().index[:down_sample])
        analyze_clusters(df=df.copy(), global_dist_matrix=normalized_dist_matrix.copy(),
                         labels=labels.copy(), n_clusters=down_sample,
                         save_path=os.path.join(save_path,f'down_sample_{down_sample}'),
                         global_onehot_dist_matrix=onehot_dist_matrix.copy(),
                         cluster_indexes=cluster_idx)
    print('doing greedy down sample')
    if greedy_down_sample is not None:
        # down sample by looking for representative sequences which
        # are most diverse from already chosen sequences
        assert threshold is not None ,'threshold must be a positive integer '
        greedy_sequence_choice_algorithm(representative_variants=representative_sequences,
                                         global_dist_matrix=normalized_dist_matrix,
                                         df=df,
                                         save_path=os.path.join(save_path,f'greedy_down_sample_{greedy_down_sample}_threshold_{threshold}'),
                                         labels=labels,
                                         threshold=threshold,
                                         greedy_down_sample=greedy_down_sample)


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


def kmeans_clustering_onehot(top_sequences_path,
                             save_path,
                             onehot_matrix_path,
                             n_clusters=5,
                             random_state=0,
                             down_sample=None,
                             percentage2include=0.15,
                             threshold =100,
                             greedy_down_sample=None
                             ):
    '''

    :param top_sequences_path:
    :param save_path:
    :param onehot_matrix_path:
    :param n_clusters:
    :param random_state:
    :param down_sample:
    :param percentage2include:
    :param threshold:
    :param greedy_down_sample:
    :return:
    '''

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


    print('starting greedy down sample')
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


    print('starting standard down sample')
    if down_sample is not None:
        cluster_idx=np.array(pd.Series(labels).value_counts().index[:down_sample])
        analyze_clusters(df=df.copy(), global_dist_matrix=global_dist_matrix,
                         labels=labels.copy(), n_clusters=down_sample,
                         save_path=os.path.join(save_path, f'down_sample_{down_sample}'),
                         global_onehot_dist_matrix=onehot_dist_matrix.copy(),
                         cluster_indexes=cluster_idx)

def analyze_clusters(df,
                     labels,
                     n_clusters,
                     save_path,
                     X=None,
                     global_dist_matrix=None,
                     global_onehot_dist_matrix=None,
                     cluster_indexes=None,
                     percentage2include=0.15):
    '''

    :param df:
    :param labels:
    :param n_clusters:
    :param save_path:
    :param X:
    :param global_dist_matrix:
    :param global_onehot_dist_matrix:
    :param cluster_indexes:
    :param percentage2include:
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


def greedy_sequence_choice_algorithm(representative_variants,
                                    global_dist_matrix,
                                     df,
                                     greedy_down_sample,
                                     save_path,
                                     labels,
                                     threshold):
    '''
    greedy algorithm to choose the
    :param representative_variants: representative variants of each cluster
    :param global_dist_matrix: np.array, shape N x N,  the entire distance numpy distance matrix
                                    from all the sequences in the dataframe (df)
    :param df: pd.DataFrame() of all the top sequences that come from all the simulated annealing runs
    :param greedy_down_sample: int, number of sequences to choose , must be less than or equal to
                                    the number of unique integers in labels array
    :param save_path: the path to save the chosen sequences
    :param labels: label each data point in df as a numpy.array
    :param threshold: the minimum number of sequences to include a cluster in a final result
    :return: save chosen sequences along with position and Amino Acid distribution to save_path directory
    '''

    print('\n ===== starting greedy approach =====')

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





def run_cluster_precomputed():
    run_params = ['rejection_10','only_train_10']
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


                clustering_object = AgglomerativeClustering(n_clusters=n_cluster,
                                                            metric='precomputed',
                                                            linkage=linkage,
                                                            compute_distances=True)

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


def preprocessing_for_clustering(run_name):
    '''
    how to preprocess the run result to get run the clustering algorithm to select the sequences
    of interest.
    1. computes distance matrix for the blosum62 distances
    2. does analysis on those distances saving the average distance to every other point histogram
    3. calculates the distance matrix in one hot format
    :param run_name: name of the run as shown by  f'3d_{run_name}_mutant_10k_run'
    :return:
    '''

    start  = time.time()
    compute_distances(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                                  'analysis','top_sequences.csv'),
                                          save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                                               'analysis'))
    second = time.time()
    print( 'compute distancres time: ',second-  start)

    dist_matrix_analysis(matrix_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                              'analysis','vectorized_blosum_similarity_distance_matrix.npy'), save_path=
                                        os.path.join('results', f'3d_{run_name}_mutant_10k_run','analysis'))

    third  = time.time()
    print('dist matrix analysis: ',third-second)
    #     ## 1)

    calculated_dist_matrix_onehot(top_sequences_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                          'analysis','top_sequences.csv'),
                                   save_path=os.path.join('results', f'3d_{run_name}_mutant_10k_run',
                                                     'analysis'))

    print('onehot dist matrix time: ',time.time()-third)



def combine_into_final_results():
    runs= ['only_train_5','rejection_5','only_train_10','rejection_10']
    run_dfs= []
    for run in runs:
        run_df = pd.read_csv(os.path.join('results', f'3d_{run}_mutant_10k_run',
                                        'analysis','agglomerate_complete','blosum62_clusters_20',
                                          'greedy_down_sample_5_threshold_100',
                                          'chosen_sequences.csv'))



        run_df['run_name']=run
        run_df['best_sequence']=run_df['best_mutant'].apply(lambda x:variant2sequence(x,av_gfp_WT))
        run_dfs.append(run_df)

    df  = pd.concat(run_dfs,ignore_index=True)


    df.to_csv(os.path.join('results','final_chosen_sequences','final_sequences_my_version.csv'))



    unique_df  = df[['run_name','unique_mutations','unique_residues']].set_index('run_name')

    unique_df = unique_df[::5]

    unique_df.plot.bar()
    plt.tight_layout()
    plt.savefig(os.path.join('results','final_chosen_sequences','final_sequences_count.png'))


    df= df[['run_name', 'best_sequence', 'best_mutant', 'best_fitness']]

    # best scoring variant
    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    data_df = pd.read_csv(os.path.join('data', 'avgfp.csv'))
    train_df = data_df.iloc[train_idx]


    best_train=  train_df.iloc[train_df['score'].argmax()]


    best_train_df =  pd.DataFrame(columns=['run_name', 'best_sequence', 'best_mutant', 'best_fitness'],
                                  data=[['best_train',variant2sequence(best_train['variant'],av_gfp_WT),best_train['variant'],best_train['score']],
                                         ['wildtype',av_gfp_WT,None,0]])

    final_df = pd.concat([df,best_train_df],ignore_index=True)


    for row in final_df.itertuples():
        assert len(row.best_sequence)==237
    final_df.to_csv(os.path.join('results','final_chosen_sequences','final_sequences.csv'))



if __name__ == '__main__':
    # combine_into_final_results()
    # run_name  = 'rejection_10'
    # preprocessing_for_clustering(run_name=run_name)
    # run_kmeans_onehot()
    run_cluster_precomputed()




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