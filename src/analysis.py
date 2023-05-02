
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
from utils import onehot,variant2sequence
import torchinfo
import torchextractor as tx
import metl
import torch



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



def preprocess_10k_run():
    D=[]
    for i in np.arange(10000):
        df= pd.read_csv(os.path.join('results','10k_full_base_run',f"{i}.csv"))
        D.append(df[df.index == 9999])
    df= pd.concat(D)
    df.to_csv(os.path.join('results','10k_full_base_run','top_sequences.csv'))
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
def kmeans_clustering_onehot(dir='10k_full_base_run'):
    df = pd.read_csv(os.path.join('results',dir, 'top_sequences.csv'))

    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())

    X = np.vstack(df['onehot'].to_numpy())

    kmeans= KMeans(n_clusters=5,random_state=0).fit(X)
    #todo: start clustering these mutations, maybe get the intertia,
    # and sizes of each cluster
    arr = kmeans.labels_
    for nb in np.arange(5):
        print(sum(arr[arr==nb]))


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



if __name__ == '__main__':
    # random_turnup()
    for nb in [5,10,15]:
        threeD_vs_oneD(nb)
    # extrapolation_1d_model_sim_anneal()
    # kmeans_clustering_onehot()
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