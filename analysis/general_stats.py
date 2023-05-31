

import numpy as np
import pandas as pd
from src.utils import AAs,av_gfp_WT
import os
import matplotlib.pyplot as plt

AA_list = np.array(list(AAs))
def amino_acid_distribution(variants,save_dir=None,suffix=None):
    '''
    look at the amino acid distribution for a given set of variants, given as an iterable object
    i.e. ['R4T','R5T,T6Y']
    :param variants: input variants, can contain many mutants
    :param save_dir: directory to save a figure (if None given no figure will save)
    :param suffix: trailing character on the file name  (aa_dist_{suffix}.png)
    :return: a pandas series with the count for each of the amino acids
    '''
    train_series = pd.Series(index=list(AAs), data=np.zeros((len(AAs),)))
    for variant in variants:
        for mutant in variant.split(','):
            aa= mutant[-1]
            train_series[aa]=train_series[aa] +1



    if save_dir is not None:
        # Creating the bar plot
        fig, ax = plt.subplots(1, 1)
        ax = train_series.plot.bar(ax=ax, rot=90, alpha=0.5)
        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Occurences in Training Set')
        ax.set_title(f'Amino Acid Distribution')

        filename = 'aa_dist'
        if suffix is not None:
            filename=filename+'_'+suffix
        fig.savefig(os.path.join(save_dir,f"{filename}.png"))
    return train_series

def residue_distribution(variants,res_range=None,save_dir=None,suffix=None):
    '''
        look at the residue distribution for a given set of variants of AVGFP ONLY, given as an iterable object
        i.e. ['R4T','R5T,T6Y']
        **
        :param variants: input variants, can contain many mutants
        :param res_range: residue range (default is the range of residues for AV_GFP domain)
        :param save_dir: directory to save a figure of histogram of residues (if None given no figure will save)
        :param suffix: trailing character on the file name  (res_dist_{suffix}.png)
        :return: a pandas series with the count for each of the residues
        '''

    if res_range is None:
        res_range= np.arange(len(av_gfp_WT))


    train_series = pd.Series(index = res_range, data=np.zeros((len(res_range),)))
    for variant in variants:
        for mutant in variant.split(','):
            idx =int(mutant[1:-1])
            train_series[idx] = train_series[idx] + 1

    if save_dir is not None:
        # Creating the bar plot
        fig, ax = plt.subplots(1, 1)
        ax = train_series.plot.hist(ax=ax,alpha=0.5,bins=len(res_range)/3)
        ax.set_xlabel('Residue')
        ax.set_ylabel('Occurences in Training Set')
        ax.set_title(f'Residue Distribution')

        filename = 'res_dist'
        if suffix is not None:
            filename = filename + '_' + suffix
        fig.savefig(os.path.join(save_dir, f"{filename}.png"))

    return train_series



def order_mutant(mutant):

    position = int(mutant[1:-1])
    aa_location = np.argmax(AA_list==mutant[-1]) / 100
    return position+ aa_location

def sort_variant(variant):
    return ','.join(sorted(variant.split(','),key=lambda x:order_mutant(x)))
def unique_sequences(variants,other_unique=2,save_dir=None,suffix=None):
    '''
    finds the unique sequences then does value_counts() and plots the distribution or a bar plot
    depending on the number of found unique sequences.
    :param variants: input variants, can contain many mutants per variant
    :param other_percentage: the percentage in which to group all variants into the same bin
    :param save_dir:  directory to save a figure of pie chart of unique sequences (if None given no figure will save)
    :param suffix: trailing character on the file name  (unique_variants_dist_{suffix}.png)
    :return: the value counts data without filter from the clustering from the other percentage
    '''

    # my unique counter is completely incorrect.
    variants = variants.apply(lambda x:sort_variant(x))
    data= variants.value_counts()


    # how many variants share a 80% of mutations?

    # bc
    # ba
    # bq
    # bc

    # Calculate the percentage values


    # Identify sectors with a percentage lower than 25%
    lower_than_other = data[data < other_unique]

    # Sum the percentages of lower sectors into a single 'other' sector
    other_percentage = lower_than_other.sum()

    # Remove the lower sectors from the percentages Series
    data.drop(lower_than_other.index, inplace=True)
    # Add the 'other' sector with the summed percentage
    data[f'other (less than < {other_unique} \n' \
         f'duplicates per seq)'] = other_percentage

    if save_dir is not None:
        fig,ax= plt.subplots(figsize=(12, 8))
        # Plot the pie chart
        ax=data.plot(kind='pie', autopct='%1.1f%%',ax=ax)
        ax.set_title('Percentage of Unique Variants')
        ax.set_axis_off()
        # Display the chart
        filename = 'unique_variants_dist'
        if suffix is not None:
            filename = filename + '_' + suffix

        fig.tight_layout()
        fig.savefig(os.path.join(save_dir, f"{filename}.png"))
    return data


if __name__ == '__main__':
    suffix= 'only_train_5'
    df= pd.read_csv(os.path.join('results',f'3d_{suffix}_mutant_10k_run','analysis','top_sequences.csv'))
    unique_sequences(df['best_mutant'],save_dir=os.path.join('results',f'3d_{suffix}_mutant_10k_run','analysis'),
                     other_unique=100)




# def aa_quandrant_sequences(varaints,save_dir=None,suffix=None):
#     filtered_stats = pd.DataFrame(index=list(AAs), columns=['I', 'II', 'III'])
#     quadrant_name = 'I'
#     for start, end in tqdm([[-np.inf, -1.5], [-1.5, 0], [0, np.inf]]):
#         filtered_df = pd.read_csv(os.path.join('data', 'distributions', f'train_{start}_{end}.csv'))
#         filtered_mutants = ",".join(filtered_df['variant']).split(',')
#         filtered_series = pd.Series(index=list(AAs), data=np.zeros((len(AAs),)))
#         for variant in filtered_mutants:
#             aa = variant[-1]
#             filtered_series[aa] = filtered_series[aa] + 1
#     pass
