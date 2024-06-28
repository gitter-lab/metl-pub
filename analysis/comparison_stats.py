import os

import numpy as np
import pandas as pd
from general_stats import amino_acid_distribution,residue_distribution
from src.utils import AAs,av_gfp_WT
import seaborn as sns
import matplotlib.pyplot as plt


def score_distribution_comparison(scores , labels , save_dir=None):
    fig,ax= plt.subplots(1,1)
    for label,score in zip(labels,scores):
        ax=score.hist(ax=ax,alpha=0.5,label=label,bins=100)

    ax.set_title('Histogram of DMS score for \n '
                 f'{labels} datasets')
    ax.set_ylabel('Count')
    ax.set_xlabel('Score')
    ax.legend()
    fig.savefig(os.path.join(save_dir,f'score_dist_{labels}.png'),dpi=300,
                    bbox_inches='tight')

def aa_distribution_comparison(variants_list,labels,save_dir=None):
    '''
    plots the comparison in distribution of the variants from different runs.
    :param variants_list:
    :param labels:
    :param save_dir:
    :return:
    '''

    df = pd.DataFrame(index=list(AAs),columns=labels)
    for label,variants in zip(labels,variants_list):
        aa_dist=  amino_acid_distribution(variants)
        #guaranteed no null values for each AA in the distribution
        df[label] = aa_dist

    if save_dir is not None:
        fig, ax = plt.subplots(1, 1)
        ax = df.plot.bar(ax=ax, rot=90, alpha=0.5)
        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Occurences')
        ax.set_title(f'Amino Acid Distribution for\n'
                     f'{labels} datasets')
        # Displaying the bar plot
        ax.legend()
        fig.savefig(os.path.join(save_dir,f"aa_dist_{labels}.png"),
                    dpi=300,
                    bbox_inches='tight')

# def res_distribution_heatmap(variants_list, labels, save_dir=None):
#     # Create a DataFrame to store the counts
#     df = pd.DataFrame(index = np.arange(len(av_gfp_WT)))
#
#     for label, variants in zip(labels, variants_list):
#         df[label] = residue_distribution(variants)
#
#     # Transpose the DataFrame for better visualization
#     # res_counts = res_counts.T
#
#     # Plotting the heatmap
#     plt.figure(figsize=(15, 8))
#     sns.heatmap(df, annot=True, linewidths=.5)
#
#     plt.title('Residue Distribution Heatmap for Different Datasets', fontsize=16)
#     plt.xlabel('Residue', fontsize=12)
#     plt.ylabel('Dataset', fontsize=12)
#
#     if save_dir is not None:
#         os.makedirs(save_dir, exist_ok=True)
#         plt.savefig(os.path.join(save_dir, f'residue_distribution_heatmap_{labels}.png'), bbox_inches='tight')

def custom_annot(value):
    if pd.notna(value) and value != 0:  # Check if not NaN and not equal to 0
        return str(value)
    return ''


def res_distribution_heatmap(variants_list, labels, label_colors, save_dir=None,fn_name=None):
    # Create a DataFrame to store the counts
    df = pd.DataFrame(index=np.arange(len(av_gfp_WT)))

    for label, variants in zip(labels, variants_list):
        df[label] = residue_distribution(variants)

    df_labels=[]
    df_markers =[]
    for i in np.arange(len(av_gfp_WT)):
        inner_labels=[]
        inner_markers=[]
        for j in np.arange(len(labels)):
            if df.iloc[i,j]==0:
                inner_labels.append('')
                inner_markers.append(0)
            else:
                inner_markers.append(j+1)
                inner_labels.append(df.iloc[i,j])
        df_labels.append(inner_labels)
        df_markers.append(inner_markers)
    df_labels = np.array(df_labels).T
    df_markers =np.array(df_markers).T

    labels_str = df_labels.astype(str)
    cmap = sns.color_palette(label_colors, as_cmap=True)

    df_final = pd.DataFrame(columns=np.arange(len(av_gfp_WT)),index=labels,data=df_markers)
    sns.heatmap(df_final, annot=labels_str, fmt='',
                linewidths=.5, cbar=False, cmap=cmap,annot_kws={"size": 3,"color": "black"})

    # plt.title('', fontsize=16)
    plt.xlabel('Position')
    plt.ylabel('Mutant')

    x_ticks = range(0, len(df_final.columns), 5)
    x_tick_labels = df_final.columns[x_ticks]
    plt.xticks(x_ticks, x_tick_labels)

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        if fn_name is not None:
            pass
        else:

            fn_name=f'res_dist_heatmap_{labels}.pdf'
        plt.savefig(os.path.join(save_dir,fn_name), bbox_inches='tight')




def res_distribution_heatmap(variants_list, labels, label_colors, save_dir=None,fn_name=None,pos_num=None):
    # Create a DataFrame to store the counts
    df = pd.DataFrame(index=np.arange(len(av_gfp_WT)))

    for label, variants in zip(labels, variants_list):
        df[label] = residue_distribution(variants)

    df_labels=[]
    df_markers =[]
    for i in np.arange(len(av_gfp_WT)):
        inner_labels=[]
        inner_markers=[]
        for j in np.arange(len(labels)):
            if df.iloc[i,j]==0:
                inner_labels.append('')
                inner_markers.append(0)
            else:
                inner_markers.append(j+1)
                if pos_num is None:
                    inner_labels.append(df.iloc[i,j])
                else:
                    inner_labels.append(pos_num[j][i])
        df_labels.append(inner_labels)
        df_markers.append(inner_markers)
    df_labels = np.array(df_labels).T
    df_markers =np.array(df_markers).T

    labels_str = df_labels.astype(str)
    cmap = sns.color_palette(label_colors, as_cmap=True)

    df_final = pd.DataFrame(columns=np.arange(len(av_gfp_WT)),index=labels,data=df_markers)
    sns.heatmap(df_final, annot=labels_str, fmt='',
                linewidths=.5, cbar=False, cmap=cmap,annot_kws={"size": 3,"color": "black"})

    # plt.title('', fontsize=16)
    if pos_num is None:
        plt.ylabel('Mutant')
    else:
        plt.ylabel('Bin')
    plt.xlabel('Position')


    x_ticks = range(0, len(df_final.columns), 5)
    x_tick_labels = df_final.columns[x_ticks]
    plt.xticks(x_ticks, x_tick_labels)

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        if fn_name is not None:
            pass
        else:

            fn_name=f'res_dist_heatmap_{labels}.pdf'
        plt.savefig(os.path.join(save_dir,fn_name), bbox_inches='tight')

    # plt.show()

def res_distribution_comparison(variants_list,labels,save_dir=None):
    # this function would take in lists of variants
    fig, ax = plt.subplots(figsize=(25, 6))
    res_dist  = pd.DataFrame()
    for label, variants in zip(labels, variants_list):
        res_dist[label]= residue_distribution(variants)

    res_dist_no_labels = res_dist.copy()
    res_dist_no_labels[res_dist_no_labels != 0] = 1
    ax = res_dist_no_labels.plot.bar(ax=ax,stacked=True, alpha=0.5)
    ax.set_title('Histogram of Residues for \n '
                 f'{labels} datasets')
    ax.set_ylabel('Count')
    ax.set_xlabel('Residue')
    ax.legend()
    # Adjust x-axis tick positions and labels
    x_ticks = range(0, len(res_dist.index), 5)
    x_tick_labels = res_dist.index[x_ticks]
    ax.set_xticks(x_ticks, x_tick_labels)

    for bar,label2patch in zip(ax.patches,res_dist.values.T.flatten()):
        if label2patch!=0:
            ax.text(
                # Put the text in the middle of each bar. get_x returns the start
                # so we add half the width to get to the middle.
                bar.get_x() + bar.get_width() / 2,
                # Vertically, add the height of the bar to the start of the bar,
                # along with the offset.
                bar.get_height()/2 + bar.get_y(),
                # This is actual value we'll show.
                label2patch,
                # Center the labels and style them a bit.
                ha='center',
                weight='bold',
                size=8
            )




    fig.savefig(os.path.join(save_dir, f'res_dist_{labels}.png'), dpi=300,
                bbox_inches='tight')

def unique_distribution_comparison(variants_list,labels,save_dir=None):
    fig, ax = plt.subplots(1, 1)
    for label, variants in zip(labels, variants_list):
        ax = variants.value_counts().hist(ax=ax, alpha=0.5, label=label, bins=100)
    ax.set_title('Histogram of Unique Sequences for \n '
                 f'{labels} datasets')
    ax.set_ylabel('Count')
    ax.set_xlabel('# of Unique Sequences')
    ax.legend()
    fig.savefig(os.path.join(save_dir, f'unique_dist_{labels}.png'), dpi=300,
                bbox_inches='tight')

    # this function would take in lists of variants
    pass
def pairwise_overlap(variants_list,labels,save_dir):
    # takes in lists of variants calculates pairwise overlap
    # between all the possible sequences given.
    df= pd.DataFrame(index=np.arange(len(labels)),columns=np.arange(len(labels)))
    data=[]
    for i,label_i in enumerate(labels):
        row=[]
        for j,label_j in enumerate(labels):
            counts1 = variants_list[i].value_counts()
            counts2 = variants_list[j].value_counts()
            common_values = counts1.index.intersection(counts2.index)
            common_count = len(common_values)
            # row then column
            row.append(float(common_count))
        data.append(row)





    fig,ax =plt.subplots()
    ax = sns.heatmap(pd.DataFrame(data=data,index=labels,columns=labels),cmap='coolwarm', linewidths=0.5, linecolor='white', annot=True, fmt='.1f', cbar=True, ax=ax)
    ax.set_title('Pairwise Overlap', fontsize=16)  # Add a title to the heatmap
    # ax.set_xlabel('X-axis Label', fontsize=12)  # Add a label to the x-axis
    # ax.set_ylabel('Y-axis Label', fontsize=12)  # Add a label to the y-axis

    # Adjust the tick labels font size
    # ax.tick_params(axis='both', which='both', labelsize=10)

    # Save or show the heatmap
    fig.tight_layout()  # Adjust the layout for better visualization
    fig.savefig(os.path.join(save_dir,f'overlap_{labels}.png'), dpi=300) # Save the heatmap as an image file




# look at the earth mover distance between the res_distribution and the aa_distribution (or KL divergence)
# or I can look at these other metrics as well, this shouldn't be too hard!!
# then I can see if I get reproducibility!

if __name__ == '__main__':
    nb_mutations= 10
    labels= [f'{nb_mutations}',f'rejection_{nb_mutations}',f'only_train_{nb_mutations}']
    variants_list = []
    scores = []
    for label in labels:
        df=pd.read_csv(os.path.join('results',f'3d_{label}_mutant_10k_run','analysis','top_sequences.csv'))
        variants_list.append(df['best_mutant'])
        scores.append(df['best_fitness'])

    # score_distribution_comparison(scores,labels,save_dir=os.path.join('results','comparisons'))
    # aa_distribution_comparison(variants_list,labels,save_dir=os.path.join('results','comparisons'))
    # unique_distribution_comparison(variants_list,labels,save_dir=os.path.join('results','comparisons'))

    label_colors = sns.color_palette("coolwarm", n_colors=3, as_cmap=False)


    res_distribution_heatmap(variants_list,labels,label_colors)
    # res_distribution_comparison(variants_list,labels, save_dir=os.path.join('results','comparisons'))
    # such a bad test to write because their is literally no overlap at all between any of the
    # situtations, but that may not be true with clustering so I should figure out that bug.
    # pairwise_overlap(variants_list,labels,save_dir=os.path.join('results','comparisons'))








