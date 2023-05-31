import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from src.utils import AAs,av_gfp_WT
from unit_tests.preprocessing_gfp_run import preprocess_train_variants_for_sim_anneal
import seaborn as sns
from tqdm import tqdm
from comparison_stats import res_distribution_comparison
def analyze_training_set():
    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    df = pd.read_csv(os.path.join('data', 'avgfp.csv'))
    train_df = df.iloc[train_idx]

    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Plot the histogram
    n_bins = 25
    color = '#b7e4cf'
    edgecolor = 'white'
    linewidth = 1.2
    ax.hist(train_df['score'], bins=n_bins, color=color, edgecolor=edgecolor, linewidth=linewidth)

    # Add labels and title
    ax.set_xlabel('Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title('Histogram of Scores in Training Set', fontsize=16, fontweight='bold')

    # Customize tick labels and grid
    ax.tick_params(axis='both', labelsize=10)
    ax.grid(axis='y', alpha=0.3)

    # Save the plot
    fig.savefig(os.path.join('data', 'distributions', 'train64_distribution.png'), dpi=300, bbox_inches='tight')
def analyze_entire_dataset():
    df = pd.read_csv(os.path.join('data', 'avgfp.csv'))

    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Plot the histogram
    n_bins = 100
    color = '#b7e4cf'
    edgecolor = 'white'
    linewidth = 1.2
    ax.hist(df['score'], bins=n_bins, color=color, edgecolor=edgecolor, linewidth=linewidth)

    # Add labels and title
    ax.set_xlabel('Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title('Histogram of Scores in Entire DataSet', fontsize=16, fontweight='bold')

    # Customize tick labels and grid
    ax.tick_params(axis='both', labelsize=10)
    ax.grid(axis='y', alpha=0.3)

    # Save the plot
    fig.savefig(os.path.join('data', 'distributions', 'entire_dataset_distribution.png'), dpi=300, bbox_inches='tight')



def amino_acids_in_training_data(normalized= False):
    # look at frequency of aa's in training set, look at frequency of amino acids in specific quadrants as well
    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    df = pd.read_csv(os.path.join('data', 'avgfp.csv'))

    train_df = df.iloc[train_idx]

    training_variants = ",".join(train_df['variant']).split(',')


    train_series = pd.Series(index=list(AAs),data=np.zeros((len(AAs),)))

    for variant in training_variants:
        aa = variant[-1]
        train_series[aa] = train_series[aa] + 1

    # Creating the bar plot
    fig,ax  = plt.subplots(1,1)
    ax= train_series.plot.bar(ax=ax,rot=90,alpha=0.5)
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Occurences in Training Set')
    ax.set_title('Training Set Avgfp Amino Acid Distribution \n'
                 '(NOT unique variants)')

    # Rotating the x-axis labels for better readability

    # Displaying the bar plot
    fig.savefig(os.path.join('data', 'distributions', 'aa_train_distribution.png'), dpi=300, bbox_inches='tight')

    # Now lets look at the distribution for each quadrant


    filtered_stats =pd.DataFrame(index = list(AAs),columns=['I','II','III'])

    quadrant_name =  'I'

    for start, end in tqdm([[-np.inf, -1.5], [-1.5, 0], [0, np.inf]]):
        filtered_df= pd.read_csv(os.path.join('data', 'distributions', f'train_{start}_{end}.csv'))
        filtered_mutants = ",".join(filtered_df['variant']).split(',')
        filtered_series = pd.Series(index=list(AAs),data=np.zeros((len(AAs),)))
        for variant in filtered_mutants:
            aa = variant[-1]
            filtered_series[aa] = filtered_series[aa] + 1


        if normalized:
            filtered_series  = filtered_series / len(filtered_mutants)
        # Creating the bar plot
        fig, ax = plt.subplots(1, 1)
        ax = filtered_series.plot.bar(ax=ax, rot=90,alpha=0.5)
        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Occurences in Training Set')
        ax.set_title(f'Training Set Avgfp Amino Acid Distribution\n'
                     f' from [{start:0.2f},{end:0.2f}]')
        # Displaying the bar plot
        fig.savefig(os.path.join('data', 'distributions', f'aa_train_distribution_normalized_{normalized}_{start}_{end}.png'), dpi=300, bbox_inches='tight')


        filtered_stats[quadrant_name] = filtered_series
        # trick to update quadrant name
        quadrant_name= quadrant_name+'I'

    fig, ax = plt.subplots(1, 1)
    ax = filtered_stats.plot.bar(ax=ax, rot=90, alpha=0.5)
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Occurences in Training Set')
    ax.set_title(f'Training Set Avgfp Amino Acid Distribution')
    # Displaying the bar plot
    ax.legend(title='Quadrant')
    fig.savefig(os.path.join('data', 'distributions', f'aa_train_distribution_normalized_{normalized}_quadrants.png'), dpi=300,
                bbox_inches='tight')



def quadrant_res_comparison():
    labels = ['I','II','III']
    variant_list= []
    for start, end in tqdm([[-np.inf, -1.5], [-1.5, 0], [0, np.inf]]):
        filtered_df= pd.read_csv(os.path.join('data', 'distributions', f'train_{start}_{end}.csv'))
        variant_list.append( filtered_df['variant'])
    res_distribution_comparison(variant_list,labels,save_dir=os.path.join('data', 'distributions'))
def number_distribution_in_training_set():
    quadrant_name='I'
    filtered_stats = pd.DataFrame(columns=['I', 'II', 'III'])
    tot=[]
    for start, end in tqdm([[-np.inf, -1.5], [-1.5, 0], [0, np.inf]]):
        filtered_df= pd.read_csv(os.path.join('data', 'distributions', f'train_{start}_{end}.csv'))
        filtered_series = pd.Series(index = np.arange(10),data=np.zeros((10,)))
        for element in filtered_df['num_mutations']:
            filtered_series[element] = filtered_series[element]+1
        filtered_stats[quadrant_name]=filtered_series
        quadrant_name = quadrant_name+ 'I'
        tot.append(len(filtered_df))


    fig, ax = plt.subplots(1, 1)
    ax = filtered_stats.plot.bar(ax=ax, rot=90, alpha=0.5)
    ax.set_xlabel('Number of Mutations')
    ax.set_ylabel('Occurences in Training Set')
    ax.set_title(f'Training Set Avgfp Number of Mutations per Variant Distribution\n'
                 f'I:{tot[0]}\n'
                 f'II:{tot[1]}\n'
                 f'III:{tot[2]}\n')
    # Displaying the bar plot
    ax.legend(title='Quadrant')
    fig.savefig(os.path.join('data', 'distributions', f'num_mutations_train_distribution_quadrants.png'), dpi=300,
                bbox_inches='tight')

def make_heat_maps_for_quadrants():
    aa_list=[]
    for aa in AAs:
        aa_list.append(aa)

    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    avgfp_df = pd.read_csv(os.path.join('data', 'avgfp.csv'))
    train_df = avgfp_df.iloc[train_idx]

    for start,end in [[-np.inf,-1.5],[-1.5,0],[0,np.inf]]:
        # populate pandas dataframe with these splits
        df = pd.DataFrame(columns=np.arange(len(av_gfp_WT)),
                          data = np.zeros((len(aa_list),len(av_gfp_WT))),
                          index=aa_list)
        filtered_bool= (train_df['score'] > start) &(train_df['score'] < end)
        filtered= train_df[filtered_bool]
        filtered.to_csv(os.path.join('data','distributions',f'train_{start}_{end}.csv'))
        for variant in filtered['variant']:
            for v in variant.split(','):
                df.loc[v[-1], int(v[1:-1])]= df.loc[v[-1], int(v[1:-1])]+ 1

        fig,ax= plt.subplots()
        ax=sns.heatmap(df, cmap='coolwarm',ax=ax)
        ax.set_title(f'Start : {start}, End: {end}, nb_variants:{len(filtered)}')
        fig.savefig(os.path.join('data','distributions',f'train_{start}_{end}.png'))




if __name__ == '__main__':
    # analyze_training_set()
    # amino_acids_in_training_data()
    number_distribution_in_training_set()

    # quadrant_res_comparison()
    # analyze_entire_dataset()

    # make_heat_maps_for_quadrants()