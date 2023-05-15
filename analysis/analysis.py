
### this script will run k-means from scikit-learn,
### also maybe where we can visualize where the sequences are
## if we run PCA! using embeddings for METL-L

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import seaborn as sns
import matplotlib.cm as cm
from src.utils import onehot,variant2sequence
from unit_tests.preprocessing_gfp_run import preprocess_train_variants_for_sim_anneal
import torchinfo
import torchextractor as tx
import metl
import torch
from tqdm import tqdm
import seaborn as sns
from scipy.spatial.distance import cdist

def random_turnup():
    '''
    Predictions for randomly generated variants.
    For increasing numbers of mutations (n=1, 2, 3, …),
    predict scores for X randomly generated variants.
    Plot the resulting prediction distributions and means.
    I would expect a downward trend in mean predicted score as we increase n.

    Progress-- 03.31.23
    Just submitted a run on a random sampler for 1-51 mutations.
    Going to take the results of this and plot them as the average prediction maybe?
    Or maybe increments of 10 as a violin plot. Or both.
    :return:
    '''

    #make the violin plots to save
    df=pd.DataFrame()
    S, F = [], []
    for uuid in np.arange(start=1,stop=50,step=5):
        temp_df = pd.read_csv(os.path.join('results','random_turnup', f"{uuid}.csv"))
        length = len(temp_df)
        F += list(temp_df['current_fitness'])
        S += [uuid] * length
    df['current_fitness'] = F
    df['nb_mutations'] = S
    fig, ax = plt.subplots(1, 1)
    ax = sns.violinplot(data=df, x="nb_mutations", y="current_fitness", alpha=0.2, ax=ax)
    ax.set_title(f'nb mutations vs current fitness on different mutation distances')
    fig.savefig(os.path.join('results','random_turnup', f'violin_plot_current_fitness.png'))

    Mean,Std =[],[]
    for uuid in np.arange(50):
        temp_df = pd.read_csv(os.path.join('results', 'random_turnup', f"{uuid}.csv"))
        Mean.append(temp_df['current_fitness'].mean())
        Std.append(temp_df['current_fitness'].std())


    fig, ax = plt.subplots(1, 1)
    ax.errorbar(np.arange(50),Mean, yerr=Std, fmt='o', color='r',alpha=0.5)
    ax.set_xlabel('nb_mutations')
    ax.set_ylabel('mean fitness')
    ax.set_title('mean fitness vs number of mutations')

    fig.savefig(os.path.join('results','random_turnup','current_fitness_mean.png'))



def preprocess_10k_run(run_saved_directory,top_sequences_directory):
    D=[]
    for i in tqdm(np.arange(10000)):
        df= pd.read_csv(os.path.join(run_saved_directory,f"{i}.csv"))
        D.append(df[df.index == 9999])
    df= pd.concat(D)
    df.to_csv(os.path.join(top_sequences_directory,'top_sequences.csv'))
def distribution_of_best_scoring(dir='10k_full_base_run'):
    df=pd.read_csv(os.path.join('results', dir,'top_sequences.csv'))
    fig,ax= plt.subplots(1,1)
    ax=df['best_fitness'].hist(ax=ax,bins=100,alpha=0.3)
    ax.set_xlabel('best_fitness')
    ax.set_title(f'distribution best_fitness {dir}')
    ax.set_ylabel('frequency')
    fig.savefig(os.path.join('results',dir,'best_sequences_dist.png'))
def distribution_with_PCA_onehot(dir='10k_full_base_run'):
    df = pd.read_csv(os.path.join('results',dir, 'top_sequences.csv'))

    df['best_sequence']= df['best_mutant'].apply(lambda x:variant2sequence(x))
    df['onehot'] =  df['best_sequence'].apply(lambda x:onehot(x).flatten())

    X= np.vstack(df['onehot'].to_numpy())
    pca=PCA(n_components=2)
    pca.fit(X)
    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)

    X_new= pca.transform(X)
    plt.scatter(X_new[:,0],X_new[:,1],s=0.3,alpha=0.3, c=df['best_fitness'],
                   cmap=cm.bwr)
    plt.xlabel('pca1')
    plt.ylabel('pca2')
    plt.title(f'pca for {dir}, with coloring on best_sequence\n'
                 f'pca.explained_variance_ratio_ :{pca.explained_variance_ratio_}')
    plt.colorbar()
    plt.savefig(os.path.join('results','10k_full_base_run','best_sequence_pca.png'))


def kmeans_clustering_onehot(top_sequences_path,save_path,n_clusters=5,random_state=0):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())

    X = np.vstack(df['onehot'].to_numpy())

    #1 - complete the kmeans clustering
    kmeans= KMeans(n_clusters=n_clusters,random_state=random_state).fit(X)
    arr = kmeans.labels_

    stats_df = pd.DataFrame(columns=['cluster','nb_points','inertia','median point','in_training_set',
                                     'train_variants'],
                            index=np.arange(n_clusters))


    training_variants = preprocess_train_variants_for_sim_anneal()
    for nb in np.arange(n_clusters):

        filtered_df = df[arr==nb]
        filtered_X = X[arr==nb]
        dist_matrix = cdist(filtered_X, filtered_X, metric='euclidean')

        # sum up the pairwise distances
        sum_dist= np.sum(dist_matrix, axis=0)

        # agrument of the median point
        argmed = np.argmin(np.abs(sum_dist - np.median(sum_dist)))

        # get the median variant
        median_variant = filtered_df.iloc[argmed]['best_mutant']

        # find the similarity of this point to the training set.
        count = 0
        tv =[]
        for v in median_variant.split(','):
            if v in training_variants:
                tv.append(v)
                count+=1

        stats_df.iloc[nb] = [nb,len(sum_dist),np.sum(sum_dist),median_variant,count,tv]
        # Create figure and axes objects
        fig, ax = plt.subplots()

        # Set plot title and axis labels
        ax.set_title(f"cluster: {nb}, nb_points:{len(sum_dist)} , \n"
                     f"inertia (total sum) : {np.sum(sum_dist)}\n,"
                     f"median point: {median_variant},\n"
                     f"in_training_set : {count},\n"
                     f"train_variants : {tv}")

        ax.set_xlabel("pairwise distance")
        ax.set_ylabel("frequency")

        # Set number of bins for histogram
        num_bins = 200

        # Plot histogram of data with specified number of bins
        n, bins, patches = ax.hist(sum_dist, bins=num_bins, color="#6495ED", alpha=0.5)

        # Customize appearance of histogram
        for i in range(len(patches)):
            # Set edge color and line width of each bar in histogram
            patches[i].set_edgecolor("#000080")
            patches[i].set_linewidth(0.5)

        # Add vertical line for mean value
        mean = np.median(sum_dist)
        ax.axvline(mean, color="#FFA500", linestyle="--", linewidth=1.5)

        # Add legend for mean value
        ax.legend([f"Mean Value = {mean:.2f}"], loc="upper right")

        # Save plot as image file
        plt.tight_layout()
        fig.savefig(os.path.join(save_path,f'cluster_{nb}.png'), dpi=300)

    # stats_df.hist()
    # plt.savefig(os.path.join(save_path,f'hist.png'))


def threeD_vs_oneD(nb_mutation):
    '''
    this function does the parity plots for the 3d vs 1d models at different rates mutation rates.
    But these aren't the same mutations crap. Wait they should be because I used the same seed! yes,
    they 100\% should be the same. Lets check.
    :return:
    '''
    df1=pd.read_csv(os.path.join('results','3d_vs_1d',f'1d_{nb_mutation}.csv')).set_index('current_mutant')
    df3=pd.read_csv(os.path.join('results','3d_vs_1d',f'3d_{nb_mutation}.csv')).set_index('current_mutant')
    df3['1d']  = df1['current_fitness']
    df3['3d'] = df3['current_fitness']

    fig,ax = plt.subplots(1,1)
    ax=df3.plot.scatter(x='1d',y='3d',s=0.05,alpha=0.3,ax=ax)
    ax.set_ylabel('3d')
    ax.set_xlabel('1d')
    ax.set_title(f'3d vs 1d prediction 10k samples for {nb_mutation} mutants')
    fig.savefig(os.path.join('results','3d_vs_1d',f'3d_vs_1d_{nb_mutation}.png'))


def sim_anneal_extrapolation_1d_100_per_mutant():

    B,M = [],[]
    for uuid,nb_mutations in zip(np.arange(11*100),np.repeat(np.arange(start=1,stop=55,step=5),100)):
        temp_df=pd.read_csv(os.path.join('results','predict_sim_anneal_extrapolation_100_simulations_per_mutant',f'{uuid}.csv'))
        B.append(temp_df.iloc[-1]['best_fitness'])
        M.append(nb_mutations)

    df=pd.DataFrame()
    df['best_fitness']=B
    df['nb_mutations']= M

    fig, ax = plt.subplots(1, 1)
    ax = sns.violinplot(data=df, x="nb_mutations", y="best_fitness", alpha=0.2, ax=ax)
    ax.set_title('100 simulations of simulated annealing\n'
                 ' at different mutation distances')
    fig.savefig(os.path.join('results',
                             'predict_sim_anneal_extrapolation_100_simulations_per_mutant',
                             'best_over_100_sims.png'))


def simple_analysis_with_original_dataset():
    ## see what mutations were suggested etc outside original dataset along with scoring distributions
    pass


def track_trajectory_of_embeddings_in_a_dataset():
    pass

def save_embeddings(dir='10k_full_base_run'):

    model, data_encoder = metl.get_from_checkpoint(os.path.join("src", "YoQkzoLD.pt"))
    return_layers= [
        'model.flatten',
        'model.backbone.2'
    ]
    extractor = tx.Extractor(model.eval(), return_layers)
    summary = torchinfo.summary(model, depth=5, verbose=1, row_settings=["var_names"])

    wt_avgfp = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQ" \
         "HDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKN" \
         "GIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

    df = pd.read_csv(os.path.join('results', dir, 'top_sequences.csv')).head(10)
    df['intermediate_layers'] =df['best_mutant'].apply(lambda x: get_embeddings(x,wt=wt_avgfp,data_encoder=data_encoder,extractor=extractor))
    df['model_out'] = df['best_mutant'].apply(lambda x:get_score(x,wt=wt_avgfp,data_encoder=data_encoder,extractor=extractor))
    df['fully_connected_embedding'] = df['intermediate_layers'].apply(lambda x:x['model.flatten'].numpy().reshape(-1))
    df['global_pooling_embedding'] = df['intermediate_layers'].apply(lambda x: x['model.backbone.2'].numpy().reshape(-1))
    gpe = np.vstack(df['global_pooling_embedding'].to_numpy())
    np.save(os.path.join('results', dir,'global_pooling_embedding.npy'),gpe)
    fce = np.vstack(df['fully_connected_embedding'].to_numpy())
    np.save(os.path.join('results', dir, 'fully_connected_embedding.npy'), fce)
def read_embedding_csv(dir='10k_full_base_run'):
    fce2  =np.load(os.path.join('results', dir, 'global_pooling_embedding.npy'))

def get_embeddings(mutants,wt,data_encoder,extractor):
    '''
    this method only works for avgfp 3d model
    :param mutants:
    :return:
    '''
    print('stop')
    if type(mutants)==str:
        mutants=[mutants]
    encoded_seqs = data_encoder.encode_variants(wt,mutants)
    with torch.no_grad():
        model_out, intermediate_out = extractor(torch.tensor(encoded_seqs), pdb_fn="1gfl_cm.pdb")

    return intermediate_out

def get_score(mutants,wt,data_encoder,extractor):
    '''
     this method only works for avgfp 3d model
     :param mutants:
     :return:
     '''
    print('stop')
    if type(mutants) == str:
        mutants = [mutants]
    encoded_seqs = data_encoder.encode_variants(wt, mutants)
    with torch.no_grad():
        model_out, intermediate_out = extractor(torch.tensor(encoded_seqs), pdb_fn="1gfl_cm.pdb")

    return model_out[0][0]


def pca_with_melt_embeddings(dir='10k_full_base_run',embedding='global_pooling_embedding'):
    df = pd.read_csv(os.path.join('results', dir, 'top_sequences.csv'))
    X = np.load(os.path.join('results', dir,f'{embedding}.npy'))
    pca = PCA(n_components=2)
    pca.fit(X)
    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)

    X_new = pca.transform(X)
    plt.scatter(X_new[:, 0], X_new[:, 1], s=0.3, alpha=0.3, c=df['best_fitness'],
                cmap=cm.bwr)
    plt.xlabel('pca1')
    plt.ylabel('pca2')
    plt.title(f'pca for {dir}, metl ({embedding}) of best sequences\n'
              f'pca.explained_variance_ratio_ :{pca.explained_variance_ratio_}')
    plt.colorbar()
    plt.savefig(os.path.join('results', dir, f'best_sequence_pca_{embedding}.png'))

def kmeans_with_metl_embeddings(dir='10k_full_base_run',embedding='global_pooling_embedding'):
    df = pd.read_csv(os.path.join('results', dir, 'top_sequences.csv'))
    X = np.load(os.path.join('results', dir, f'{embedding}.npy'))


def extrapolation_1d_model_sim_anneal():
    '''
    keep in mind this extrapolation contains all the same seed. I don't think that matters. Since
    at different mutation distances, the sequences should diverge on the accpetance, etc.
    But its something to keep in mind for the future when looking at this analysis.

    Also keep in mind this isn't hill climbing. This is our version of sim anneal. And as we
    see in the case of the 10k run, their can be large distributions in the predicted scores.

    :return:
    '''
    F,S,B = [],[],[]
    df=pd.DataFrame()
    for uuid,nb_mutations in zip(np.arange(11),np.arange(start=1,stop=55,step=5)):
        temp_df= pd.read_csv(os.path.join('results','predict_sim_anneal_extrapolation',f'{uuid}.csv'))
        length=len(temp_df)
        F += list(temp_df['current_fitness'])
        S += [nb_mutations] * length
        B.append(temp_df.iloc[-1]['best_fitness'])



    df['fitness'] = F
    df['nb_mutations'] = S

    fig, ax = plt.subplots(1, 1)
    ax = sns.violinplot(data=df, x="nb_mutations", y="fitness", alpha=0.2, ax=ax)
    ax.set_title(f'steps in run : {length}, simulated annealing\n'
                 f'current_fitness all with same seed, 1d model ')
    fig.savefig(os.path.join('results', 'predict_sim_anneal_extrapolation','current_fitness_violin_plot_distributions.png'))


    fig, ax = plt.subplots(1, 1)
    ax.scatter(x=np.arange(start=1,stop=55,step=5),y=B,alpha=0.5)
    ax.set_xlabel('nb_mutations')
    ax.set_ylabel('best_fitness from run')
    ax.set_title(f'steps in run : {length}, simulated annealing\n'
                 f'best_fitness all with same seed, 1d model ')
    fig.savefig(os.path.join('results', 'predict_sim_anneal_extrapolation','best_fitness_trend.png'))


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

def analyze_run(uuid):
    df = pd.read_csv(os.path.join('results', f'{uuid}.csv'))
    fig, ax = plt.subplots(1, 1)

    ax = df.plot(y=['best_fitness', 'current_fitness'], ax=ax)
    ax.set_xlabel('step')
    ax.set_ylabel('fitness')
    ax.set_title('fitness vs step sim anneal')
    fig.legend()
    fig.savefig(os.path.join('results', f'{uuid}_index.png'))

    fig,ax=plt.subplots(1,1)
    ax.semilogx(df['temperature'],df['best_fitness'],label='best_fitness')
    ax.semilogx(df['temperature'],df['current_fitness'],label='current_fitness')
    ax.set_xlabel('temperature')
    ax.set_ylabel('fitness')
    ax.set_title('fitness vs temperature for sim anneal run')
    fig.legend()
    fig.savefig(os.path.join('results',f'{uuid}_temp.png'))

def violin_plots_across_runs(list_of_uuids,tags):
    df=pd.DataFrame()
    length = None
    S,F=[],[]

    for uuid,tag in zip(list_of_uuids,tags):
        temp_df= pd.read_csv(os.path.join('results',f"{uuid}.csv"))

        if length is None:
            length = len(temp_df)
        else:
            assert length == len(temp_df) ,'this function doesnt support variable length runs'

        F+=list(temp_df['current_fitness'])
        S+=[tag]*length

    df['fitness'] = F
    df['sampler'] = S

    fig,ax=plt.subplots(1,1)
    ax=sns.violinplot(data=df, x="sampler", y="fitness",alpha=0.2,ax=ax)
    ax.set_title(f'steps in run : {length}')
    fig.savefig(os.path.join('results',f'analysis_{list_of_uuids}.png'))

    print('complete!')

def get_starting_temperature(uuid,WT_fitness=0):
    df=pd.read_csv(os.path.join('results',f"{uuid}.csv"))
    print(f"lowest fitness {df['current_fitness'].min()}")

    print(f"delta E = {df['current_fitness'].min()-WT_fitness }")

def compare_plots(uuids):
    for uuid in uuids:
        df=pd.read_csv(os.path.join('results',f"{uuid}.csv"))
        print(uuid, ',max value')
        print(df['current_fitness'].max())
        print('')
def compare_to_original(uuids=None):
    df_original = pd.read_csv(os.path.join('data', 'avgfp.tsv'), sep='\t')
    fig,ax=plt.subplots(1,1)
    ax=df_original['score'].hist(alpha=0.4,bins=50,ax=ax)
    fig.show()


def best_and_current_averages(nb_mutations,run_saved_directory,save_directory):
    D=[]
    for i in tqdm(np.arange(10000)):
        sub_df= pd.read_csv(os.path.join(run_saved_directory,f"{i}.csv"))
        D.append(sub_df[['current_fitness','best_fitness']])

    df=pd.concat(D)
    df_avg = df.groupby(level=0).mean()
    df_std = df.groupby(level=0).std()
    df_avg = df_avg.rename(columns={'current_fitness':'current_avg',"best_fitness":'best_avg'})
    df_std = df_std.rename(columns={'current_fitness':'current_std','best_fitness':'best_std'})


    sns.set_palette("colorblind")

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    x = np.arange(10001)  # 10001 to include
    ax.plot(x, df_avg['current_avg'], alpha=0.5, label='Current Average')
    ax.plot(x, df_avg['best_avg'], alpha=0.5, label='Best Average')
    ax.fill_between(x, df_avg['current_avg'] - df_std['current_std'], df_avg['current_avg'] + df_std['current_std'],
                    color='grey', alpha=0.1)

    ax.fill_between(x, df_avg['best_avg'] - df_std['best_std'], df_avg['best_avg'] + df_std['best_std'],
                    color='grey', alpha=0.1)

    ax.set_xlabel('Iterations', fontsize=14)
    ax.set_ylabel('Fitness', fontsize=14)
    ax.set_title(f"{nb_mutations} Mutants Fitness Trend", fontsize=16)
    ax.grid(False)
    fig.legend(loc='center right', fontsize=12)
    fig.tight_layout()

    fig.savefig(os.path.join(save_directory, 'fitness_trend_fig.png'), dpi=300)



def count_number_of_training_variants(mutants,training_variants):
    count=0
    for mutant in mutants.split(','):
        if mutant in training_variants:
            count+=1
    return count

def mutations_in_training_data(nb_mutations,run_saved_directory,save_directory):
    training_variants= preprocess_train_variants_for_sim_anneal()

    D = []
    for i in tqdm(np.arange(10000)):
        sub_df = pd.read_csv(os.path.join(run_saved_directory, f"{i}.csv"))
        sub_df['current_train_mutants']= sub_df['current_mutant'].apply(lambda x: count_number_of_training_variants(x,training_variants))
        sub_df['best_train_mutants'] = sub_df['best_mutant'].apply(
            lambda x: count_number_of_training_variants(x, training_variants))
        D.append(sub_df[['current_train_mutants', 'best_train_mutants']])

    df=pd.concat(D)
    df_avg = df.groupby(level=0).mean()
    df_std = df.groupby(level=0).std()
    df_avg = df_avg.rename(columns={'current_train_mutants':'current_avg',"best_train_mutants":'best_avg'})
    df_std = df_std.rename(columns={'current_train_mutants':'current_std','best_train_mutants':'best_std'})

    sns.set_palette("colorblind")

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    x = np.arange(10001)  # 10001 to include
    ax.plot(x, df_avg['current_avg'], alpha=0.5, label='Current Average')
    ax.plot(x, df_avg['best_avg'], alpha=0.5, label='Best Average')
    ax.fill_between(x, df_avg['current_avg'] - df_std['current_std'], df_avg['current_avg'] + df_std['current_std'],
                    color='grey', alpha=0.1)

    ax.fill_between(x, df_avg['best_avg'] - df_std['best_std'], df_avg['best_avg'] + df_std['best_std'],
                    color='grey', alpha=0.1)

    ax.set_xlabel('Iterations', fontsize=14)
    ax.set_ylabel('Number of Mutations in Training Set', fontsize=14)
    ax.set_title(f"{nb_mutations} Mutants Similarity to Training Set Trend", fontsize=16)
    ax.grid(False)
    fig.legend(loc='center left' ,fontsize=12)
    fig.tight_layout()

    fig.savefig(os.path.join(save_directory, 'training_set_trend_fig.png'), dpi=300)

def acceptance_ratio(run_saved_directory,save_directory,nb_steps_to_average_over):
    # todo: code doesn't have an abstarction for the number of files in the directory.
    #       should be user input
    # define acceptance ratio as every 10, then 100 steps

    # loop over each dataset then calculate the acceptance ratio
    # as number of times that current_seq changes/nb_steps_to_average_over
    # put that list of acceptance ratio into a dataframe
    D = []
    for i in tqdm(np.arange(10000)):
        sub_df = pd.read_csv(os.path.join(run_saved_directory, f"{i}.csv"))
        # remove first element to give good division factors, yes
        # this will give slightly different results
        series = sub_df['current_mutant'][1:]

        # Set the increment value
        increment = nb_steps_to_average_over

        # Calculate the number of groups
        num_groups = len(series) // increment

        # Reshape the Series into a 2D array
        reshaped_series = series.values[:num_groups * increment].reshape(num_groups, increment)

        # Calculate the frequency of element changes
        change_frequency = (reshaped_series[:, :-1] != reshaped_series[:, 1:]).mean(axis=1)

        # append to global dataframe to do std and mean later
        D.append(pd.DataFrame(data=change_frequency,columns=['change_frequency']))

    df = pd.concat(D)
    # average all the acceptance ratios
    df_avg = df.groupby(level=0).mean()
    df_std = df.groupby(level=0).std()
    df_avg = df_avg.rename(columns={'change_frequency':'change_avg'})
    df_std = df_std.rename(columns={'change_frequency':'change_std'})

    # make a fancy plot out of the average acceptance ratio along with
    # error bars for the standard deviation of the acceptance ratio
    x = np.arange(len(df_avg))

    sns.set_palette("colorblind")

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    ax.plot(x, df_avg['change_avg'], alpha=0.5, label='Acceptance Ratio')
    ax.fill_between(x, df_avg['change_avg'] - df_std['change_std'],
                    df_avg['change_avg'] + df_std['change_std'],
                    color='grey', alpha=0.1)

    ax.set_xlabel(f'Acceptance Ratio every {nb_steps_to_average_over} steps', fontsize=14)
    ax.set_ylabel('Acceptance Ratio', fontsize=14)
    ax.set_title(f"Acceptance Ratio of Mutants Averaged Over {nb_steps_to_average_over} Steps\n"
                 f"for 10,000 simulations", fontsize=16)
    ax.grid(False)
    # fig.legend(loc ='center right',fontsize=12)

    fig.tight_layout()
    fig.savefig(os.path.join(save_directory, f'acceptance_ratio_{nb_steps_to_average_over}.png'), dpi=300)

if __name__ == '__main__':
    nb_mutations= 5
    #1)
    # preprocess_10k_run(run_saved_directory=os.path.join('results','3d_10_mutant_10k_run'),
    #                   top_sequences_directory=os.path.join('results','3d_10_mutant_10k_run','analysis'))

    #2) nb_mutations= 5
    # best_and_current_averages(nb_mutations=nb_mutations , run_saved_directory= os.path.join('results',f'3d_10k_run'),
    #                    save_directory=os.path.join('results',f'3d_10k_run','analysis'))

    #3)

    # mutations_in_training_data(nb_mutations=nb_mutations, run_saved_directory=os.path.join('results', f'3d_{nb_mutations}_mutant_10k_run'),
    #                     save_directory=os.path.join('results',f'3d_{nb_mutations}_mutant_10k_run','analysis'))


    #4)
    # kmeans_clustering_onehot(top_sequences_path=os.path.join('results', f'3d_{nb_mutations}_mutant_10k_run',
    #                                                          'analysis','top_sequences.csv'),
    #                          save_path=os.path.join('results', f'3d_{nb_mutations}_mutant_10k_run',
    #                                                 'analysis'))


    #5)
    acceptance_ratio(run_saved_directory=os.path.join('results',f'3d_{nb_mutations}_mutant_10k_run',),
                     save_directory=os.path.join('results',f'3d_{nb_mutations}_mutant_10k_run','analysis'),
                     nb_steps_to_average_over=10)

    # # random_turnup()
    # for nb in [5,10,15]:
    #     threeD_vs_oneD(nb)
    # extrapolation_1d_model_sim_anneal()

    # save_embeddings()
    # compare_to_original()
    # preprocess_10k_run()
    # distribution_of_best_scoring()
    # distribution_with_PCA_onehot()
    # analyze_embeddings()
    # pca_with_melt_embeddings()
    # save_embeddings()
    # sim_anneal_extrapolation_1d_100_per_mutant()
    # save_embeddings()
    # read_embedding_csv()
    # pca_with_melt_embeddings(embedding='fully_connected_embedding')
    # analyze_run('practice_run_9')
    # violin_plots_across_runs(['practice_run_7','practice_run_8'],tags=['log-no warmup','log-with warmup'])
    # compare_plots(['practice_run_5','practice_run_8','practice_run_9'])
    # 'practice_run_5','practice_run_6'],
    #                          ['linear','log' 'random','log_random_startup'])
    # get_starting_temperature('practice_run_5')


    # df_original.to_csv(os.path.join('data','avgfp.csv'))


    # print('stop')