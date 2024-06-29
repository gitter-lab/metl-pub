import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from comparison_stats import res_distribution_comparison, res_distribution_heatmap
from src.utils import  av_gfp_WT
from matplotlib.text import Text
paper_small = {
    'font.size': 6.0,
    'axes.labelsize': 6.0,
    'axes.labelweight': "bold",
    'axes.titlesize': 6.0,
    'axes.titleweight': "bold",
    'figure.titlesize': 6.0,
    'xtick.labelsize': 6.0,
    'ytick.labelsize': 6.0,
    'legend.fontsize': 6.0,
    'axes.linewidth': 1.5,
    'grid.linewidth': 1.0,
    'lines.linewidth': 1.0,
    'lines.markersize': 3.0,
    'patch.linewidth': 1.0,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.minor.width': 1.2,
    'ytick.minor.width': 1.2,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.minor.size': 1.0,
    'ytick.minor.size': 1.0,
    'legend.title_fontsize': None,
}

plt.rcParams['figure.figsize'] = [6, 4]
plt.rcParams['figure.dpi'] = 80
plt.rcParams['savefig.dpi'] = 100
plt.rcParams['figure.figsize'] = [6, 4]
plt.rcParams['figure.dpi'] = 80
plt.rcParams['savefig.dpi'] = 100
def my_plot(x, y,
            figsize=(4, 3),
            plotting_context=None,
            save_fn=None):
    '''
    practice plotting in Sam's plotting scheme.
    :param x:
    :param y:
    :param figsize:
    :param plotting_context:
    :param save_fn:
    :return:
    '''
    with sns.plotting_context(plotting_context):
        fig, ax = plt.subplots(1, figsize=figsize)

        sns.scatterplot(x=x, y=y, ax=ax)

        plt.title('Random Scatter Plot')
        plt.xlabel('Variable X')
        plt.ylabel('Variable Y')
        plt.grid(True)

        plt.tight_layout()

        tight_bbox = False
        if save_fn is not None:
            bbox_inches = "tight" if tight_bbox else None
            pad_inches = 0.025 if tight_bbox else 0.1
            plt.savefig(save_fn, dpi=300, facecolor="white", transparent=False,
                        bbox_inches=bbox_inches, pad_inches=pad_inches)

        plt.show()
        plt.close(fig)

def analyze_training_set(plotting_context):
    '''
    Analysis of the training data in histogram format.
    :param plotting_context:
    :return:
    '''
    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    df = pd.read_csv(os.path.join('data', 'avgfp.csv'))
    train_df = df.iloc[train_idx]
    label_colors = sns.color_palette("coolwarm", n_colors=3, as_cmap=False)
    with sns.plotting_context(plotting_context):
    # Set up the figure and axes

        fig, ax = plt.subplots(figsize=(6,6))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # Plot the histogram
        n_bins = 30
        color = '#b7e4cf'
        edgecolor = 'white'
        linewidth = 1.2


        data =train_df['score']
        mask_below_minus_1_5 = data < -1.5
        mask_above_9 = data > 0
        label_colors = sns.color_palette("coolwarm", n_colors=3, as_cmap=False)

        # ax.hist(, bins=n_bins, color=color, edgecolor=edgecolor, linewidth=linewidth)

        ax.hist(data[~mask_below_minus_1_5 & ~mask_above_9], bins=n_bins, alpha=0.9, color=label_colors[1],range=(-2.5,.5))
        ax.hist(data[mask_below_minus_1_5], bins=n_bins, alpha=0.9, color=label_colors[0],range=(-2.5,.5))
        ax.hist(data[mask_above_9], bins=n_bins, alpha=0.9, color=label_colors[2],range=(-2.5,.5))


    # Add labels and title
        ax.set_xlabel('Score', fontsize=12, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        # ax.set_title('Histogram of Scores in Training Set', fontsize=16, fontweight='bold')

        plt.axvline(x=-1.5, color='grey', linestyle='--', alpha=0.6, )
        plt.axvline(x=0, color='grey', linestyle='--', alpha=0.6)
        # plt.axvline(x=1.5, color='grey', linestyle='--', alpha=0.3, label='III')

        # Add labels near the top
        q=10
        plt.text(-2.,q, 'I', color=label_colors[0], alpha=0.6, ha='center', va='bottom', fontweight='bold',fontsize=24)
        plt.text(-1, q, 'II', color=label_colors[1], alpha=0.6, ha='center', va='bottom', fontweight='bold',fontsize=24)
        plt.text(.25, q, 'III', color=label_colors[2], alpha=0.6, ha='center', va='bottom', fontweight='bold',fontsize=24)

    # Customize tick labels and grid
        ax.tick_params(axis='both', labelsize=10)
        ax.grid(axis='y', alpha=0.3)

        # Save the plot
        fig.savefig(os.path.join('data', 'distributions', 'train64_distribution.pdf'), bbox_inches='tight')

def sim_anneal_extrapolation_1d_100_per_mutant(plotting_context):
    '''
    Simulated annealing extrapolation with 100 simulations at each position.
    :param plotting_context:
    :return:
    '''
    with sns.plotting_context(plotting_context):

        B,M = [],[]
        for uuid,nb_mutations in zip(np.arange(11*100),np.repeat(np.arange(start=1,stop=55,step=5),100)):
            temp_df=pd.read_csv(os.path.join('results','predict_sim_anneal_extrapolation_100_simulations_per_mutant',f'{uuid}.csv'))
            B.append(temp_df.iloc[-1]['best_fitness'])
            M.append(nb_mutations)

        df=pd.DataFrame()
        df['best_fitness']=B
        df['nb_mutations']= M

        fig, ax = plt.subplots(1,1)
        ax = sns.violinplot(data=df, x="nb_mutations", y="best_fitness", alpha=0.2, ax=ax)
        ax.set_xlabel('Number of mutations')
        ax.set_ylabel('Predicted Brightness')
        # ax.set_title('100 simulations of simulated annealing\n'
        #              ' at different mutation distances')
        fig.savefig(os.path.join('results',
                                 'predict_sim_anneal_extrapolation_100_simulations_per_mutant',
                                 'best_over_100_sims.pdf'))

def random_turnup(plotting_context):
    '''
    Predictions for randomly generated variants.
    For increasing numbers of mutations (n=1, 2, 3, …),
    predict scores for X randomly generated variants.
    Plot the resulting prediction distributions and means.
    I would expect a downward trend in mean predicted score as we increase n.
    :return:
    '''
    with sns.plotting_context(plotting_context):
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
        ax.set_xlabel('Number of mutations')
        ax.set_ylabel('Predicted Brightness')
        # ax.set_title(f'Number of Random Mutations vs Fitness')
        fig.savefig(os.path.join('results','random_turnup', f'violin_plot_current_fitness.pdf'))
def make_residue_maps_for_quadrants(plotting_context,figsize):
    '''
    Show in quadrant form where mutations are in GFP training data.
    :param plotting_context:
    :param figsize:
    :return:
    '''
    quadrant_name = 'I'
    comparison_labels  = [ ]
    for start, end in tqdm([[-np.inf, -1.5], [-1.5, 0], [0, np.inf]]):
        filtered_df= pd.read_csv(os.path.join('data', 'distributions', f'train_{start}_{end}.csv'))
        variants_list = []
        for variant in filtered_df['variant']:
            variants_list.append(pd.Series([variant]))
        comparison_labels = np.arange(len(filtered_df['variant']))+len(comparison_labels)
        with sns.plotting_context(plotting_context):
            fig, ax = plt.subplots(1, figsize=figsize)
            res_distribution_comparison(variants_list,comparison_labels,os.path.join('data', 'distributions'))

def count_residues(filter_df):
    '''

    :param filter_df:
    :return:
    '''
    var_dict={}
    for i in np.arange(len(av_gfp_WT)):
        var_dict[i]=0
    # var_dict = {:np.arange(len(av_gfp_WT))*0}
    for variant in filter_df['variant']:
        for mutant in variant.split(","):
            pos= int(mutant[1:-1])
            var_dict[pos] +=1

    return  pd.Series(var_dict)
def custom_annot(x, **kwargs):
        return f'$\\bf{{{x}}}$'

def final_sequences_total_map(plotting_context,figsize,save_dir=None,start_pos=None,end_pos=None):
    # Create a DataFrame to store the counts
    with sns.plotting_context(plotting_context):
        # plt.rc('text', usetex=True)
        # plt.rc('text.latex', preamble=r'\usepackage{lmodern}')
        fig, ax = plt.subplots(1, figsize=figsize)
        df_final_sequences = pd.read_csv('results/random_observed_unobserved/random_baseline.csv')
        df_final_sequences=df_final_sequences[:21]


        # label_colors = sns.color_palette("coolwarm", n_colors=3, as_cmap=False)
        label_colors = ['#FFFFFF','honeydew','aliceblue','papayawhip']


        if start_pos==None:
            index = np.arange(len(av_gfp_WT))
        else:
            index=np.arange(start_pos,end_pos,1)
        df = pd.DataFrame(index=index)

        df_labels = []
        df_markers = []

        colormap= {'blank':0,'wt':1,'rejection':2,'only':3}

        df_markers_inner=[]
        df_labels_inner=[]
        for i in index:
            df_markers_inner.append(colormap['blank'])
            df_labels_inner.append(av_gfp_WT[i])

        df_labels.append(df_labels_inner)
        df_markers.append(df_markers_inner)



        for j in np.arange(len(df_final_sequences.index)-1):
            inner_labels = []
            inner_markers = []


            for i in index:
                inner_labels.append('')
                inner_markers.append(colormap['blank'])

            # run_type =df_final_sequences.iloc[j]['run_name'].split('_')[0]

            for mutant in df_final_sequences.iloc[j]['variant'].split(','):
                pos,sub = int(mutant[1:-1]),mutant[-1:]

                if start_pos is not None:
                    if pos  < start_pos or pos >=end_pos:
                        pass
                    else:
                        inner_labels[pos-start_pos]=sub+'!'
                        inner_markers[pos-start_pos]=colormap['blank']
                else:
                    inner_labels[pos] =sub+'!'
                    inner_markers[pos] = colormap['blank']

            df_labels.append(inner_labels)
            df_markers.append(inner_markers)
        df_labels = np.array(df_labels)
        df_markers = np.array(df_markers)

        labels_str = df_labels.astype(str)
        df_final =pd.DataFrame(labels_str)
        # labels_str = df_final.applymap(lambda x: f"\textbf{{{x}}}")


        # formatted_labels = labels_str.applymap(lambda x: f'$\\bf{{{x}}}$')
        cmap = sns.color_palette(label_colors, as_cmap=True)

        df_final = pd.DataFrame(columns=index, index=np.arange(len(df_final_sequences.index)), data=df_markers)
        heatmap= sns.heatmap(df_final, annot=labels_str, fmt='',
                    linewidths=.5, cbar=False, cmap=cmap, annot_kws={"size": 5, "color": "black"})

        for i,item in enumerate(heatmap.axes.get_children()):
            if isinstance(item, Text):
                text= item.get_text()
                print(text)
                if len(text)>0:
                    if text[-1] == '!':
                        item.set_weight('bold')
                        item.set_text(text[:-1])





        # plt.title('', fontsize=16)

        plt.ylabel('Variant')
        # plt.legend()
        plt.xlabel('Position')


        x_ticks = np.arange(0, len(df_final.columns), 5)
        x_tick_labels = df_final.columns[x_ticks]
        x_ticks=x_ticks+0.5
        plt.xticks(x_ticks, x_tick_labels)

        y_ticks= np.arange(0,21,1)+.5
        # y_ticks_labels = ['wildtype','O5','O5','O5','O5','O5','U5','U5','U5','U5','U5',
        #                  'O10','O10', 'O10','O10','O10',
        #                  'U10','U10','U10','U10','U10']
        y_ticks_labels=['wildtype','U5','U5','U5','U5','U5','U10','U10','U10','U10','U10',
                        'O5','O5','O5','O5','O5','O10','O10', 'O10','O10','O10']

        plt.yticks(y_ticks,y_ticks_labels,rotation=0)
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            if start_pos is not None:
                fn_name= f'final_sequences_{start_pos}_{end_pos}.pdf'
            else:
                fn_name = f'final_sequences.pdf'
            # if fn_name is not None:
            #     pass
            # else:


            plt.savefig(os.path.join(save_dir,fn_name), bbox_inches='tight')
def quadrant_res_comparison(plotting_context,figsize):
    labels = ['I','II','III']
    variant_list= []
    pos_num=[]
    # maximum_num = -1
    # minimum_num = np.inf
    for start, end in tqdm([[-np.inf, -1.5], [-1.5, 0], [0, np.inf]]):
        filtered_df= pd.read_csv(os.path.join('data', 'distributions', f'train_{start}_{end}.csv'))
        filter_pos_num = count_residues(filtered_df)
        filter_pos_num[filter_pos_num==1]=''
        # if np.min(filter_pos_num)< minimum_num:
        #     minimum_num=np.min(filter_pos_num)
        # if np.max(filter_pos_num)> maximum_num:
        #     maximum_num=np.max(filter_pos_num)
        pos_num.append(filter_pos_num)
        variant_list.append( filtered_df['variant'])

    with sns.plotting_context(plotting_context):
        fig, ax = plt.subplots(1, figsize=figsize)
        label_colors = sns.color_palette("coolwarm", n_colors=3, as_cmap=False)
        label_colors=['#FFFFFF']+label_colors
        res_distribution_heatmap(variant_list,labels, label_colors,
                                 save_dir=os.path.join('data', 'distributions'),pos_num = pos_num)

def final_choosen_sequences(plotting_context,figsize):
    df = pd.read_csv('results/final_chosen_sequences/final_sequences.csv')
    fields=['only_train_5','rejection_5','only_train_10','rejection_10']
    palettes = ['coolwarm','coolwarm','coolwarm','coolwarm']
    for field,pallette in zip(fields,palettes):
        with sns.plotting_context(plotting_context):
            fig, ax = plt.subplots(1, figsize=figsize)
            label_colors = sns.color_palette(pallette, n_colors=5, as_cmap=False)
            label_colors = ['#FFFFFF'] + label_colors
            labels = np.arange(5)+1
            sub_df = df[df['run_name']==field]
            variant_list=[]
            for mutant in sub_df['best_mutant']:
                variant_list.append([mutant])

            res_distribution_heatmap(variant_list, labels, label_colors,
                                     save_dir=f'results/final_chosen_sequences',
                                     fn_name=f'res_dist_heatmap_{field}.pdf')








if __name__ == '__main__':
    # np.random.seed(1)
    # x = np.random.randn(100)
    # y = np.random.randn(100) + x
    # my_plot(x, y,plotting_context=paper_small)

    # make_residue_maps_for_quadrants(paper_small,figsize=(4,10))
    # quadrant_res_comparison(paper_small,figsize=(7.086,1))

    final_sequences_total_map(paper_small,figsize=(7.086,3),save_dir='results/random_observed_unobserved',
                              start_pos=161,end_pos=237)
    # final_choosen_sequences(paper_small,figsize=(7.086,4))
    # analyze_training_set(paper_small)
    # sim_anneal_extrapolation_1d_100_per_mutant(paper_small)
    # random_turnup(paper_small)
