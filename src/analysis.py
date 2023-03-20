
### this script will run k-means from scikit-learn,
### also maybe where we can visualize where the sequences are
## if we run PCA! using embeddings for METL-L

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
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
    compare_to_original()
    # analyze_run('practice_run_9')
    # violin_plots_across_runs(['practice_run_7','practice_run_8'],tags=['log-no warmup','log-with warmup'])
    # compare_plots(['practice_run_5','practice_run_8','practice_run_9'])
    # 'practice_run_5','practice_run_6'],
    #                          ['linear','log' 'random','log_random_startup'])
    # get_starting_temperature('practice_run_5')


    # df_original.to_csv(os.path.join('data','avgfp.csv'))


    print('stop')