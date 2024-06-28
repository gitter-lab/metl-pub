
import numpy as np
import os
import pandas as pd
from src.utils import AAs,av_gfp_WT
import random
wt= 'skgeelftgvvpilveldgdvnghkfsvsgegegdatygkltl\
kficttgklpvpwptlvttlsygvqcfsrypdhmkqhdffksampegyv\
qertiffkddgnyktraevkfegdtlvnrielkgidfkedgnilghkley\
nynshnvyimadkqkngikvnfkirhniedgsvqladhyqqntpigdgpvll\
pdnhylstqsalskdpnekrdhmvllefvtaagithgmdelyk'

def preprocess_train_variants_for_sim_anneal():
    '''
    get only the training varaints that were used in metl publication for
    experimental designs.
    :return: return only the unique 209 mutations
    '''
    # this function preprocesses the training variants for simulated annealing
    # it excludes certain
    train_idx = np.loadtxt(os.path.join('data','train_64_variants_avgfp.txt'),dtype=int)
    df= pd.read_csv(os.path.join('data','avgfp.csv'))

    train_df = df.iloc[train_idx]

    rejections = ",".join(train_df['variant'])
    all_mutants=  rejections.split(',')
    print(f"all mutants: {len(all_mutants)}")
    unique_mutants  =np.unique(all_mutants)
    print(f"unique mutants : {len(unique_mutants)}")
    for mutant in  unique_mutants:
        assert mutant in all_mutants, 'why would this not be true'
    return ','.join(unique_mutants)

def check_output_for_constraints():
    '''
    post-processing step to make sure that only
    mutants that are training mutants are iin each of the testing constraints .
    :return: True if test pasing
    '''
    rejected_mutants = preprocess_train_variants_for_sim_anneal().split(',')
    df = pd.read_csv(os.path.join('results','testing_constraints.csv'))
    current_mutants= np.unique(",".join(df['current_mutant']).split(','))
    for mutant in rejected_mutants:
        assert mutant not in current_mutants, 'a mutant should not be in the allowed current mutants'
    return True # if everything passes this test is passing

def check_constraints_all_mutants():
    '''
    Post-processing step to make sure that the below code was correct.
    ====
    note to self I went and copied this code into the sampler after the all_mutation_library
    with open(os.path.join('results','test_constraints_all_mutants.txt'),'w') as f:
        f.write(','.join(all_mutation_library))
    using that I then checked to see if any of the mutants had slipped into the solution and none had.
    :return:None
    '''
    rejected_mutants = preprocess_train_variants_for_sim_anneal().split(',')
    with open(os.path.join('results', 'test_constraints_all_mutants.txt'), 'r') as f:
        allowed_mutants = f.read().split(',')
    for mutant in rejected_mutants:
        assert mutant not in allowed_mutants , 'a not allowed mutant should not be in the allowed mutant generation list'

    assert len(rejected_mutants)+len(allowed_mutants) == 19*len(wt), '19 mutations at each residue * length is the total number of mutations'

    # okay I feel very confident now that the code is behaving as I would expect it too.
def get_all_singles(wild_type):
    '''
    get all single mutants of wildtype sequence.
    :param wild_type: str sequence
    :return:
    '''
    M = [ ]
    for i,element in enumerate(wild_type):
        for aa in AAs:
            if aa != element:
                M.append(f"{element}{i}{aa}")
    assert len(wild_type)*19== len(M)
    return M


def preprocess_only_allow_training():
    '''
    return the variants which are not in the training set.
    :return:
    '''
    training_variants = preprocess_train_variants_for_sim_anneal().split(',')
    singles= get_all_singles(wt.upper())

    exlusion_list = []
    for variant in singles:
        if variant in training_variants:
            pass
        else:
            exlusion_list.append(variant)
    assert len(exlusion_list)+209  == len(singles), 'this the hard coded truth'
    return ",".join(exlusion_list)

def print_high_bin_variants():
    '''
    save only the 64 training varaints.
    :return:
    '''
    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    df = pd.read_csv(os.path.join('data', 'avgfp.csv'))
    train_df = df.iloc[train_idx]
    train_df.drop(columns=['Unnamed: 0'], inplace=True)
    train_df.to_csv('data/train_64_variants_avgfp.csv',index=False)



def variant2seq(variants,wt):
    # must be zero based indexing in the variant
    new_sequence= list(wt)
    for mutant in variants.split(','):
        print(mutant)
        wt_aa,pos,mut_aa = mutant[0],int(mutant[1:-1]),mutant[-1]
        assert wt[pos]==wt_aa and new_sequence[pos]==wt_aa,'this must be true'
        new_sequence[pos] = mut_aa
    return ''.join(new_sequence)
def add_variant(mutants2sample,nb_mutations,Sequences,Variants,Mutations,training_variants,
                training_mutants=None):
    variant=[]
    sampled_positions=[]
    for i in range(nb_mutations):
        sampled_mutant =random.sample(mutants2sample, 1)

        sampled_positions.append(sampled_mutant[0][1:-1])
        next_round_mutants=[]
        for mutant2sample in mutants2sample:
            if mutant2sample[1:-1] not in sampled_positions:
                next_round_mutants.append(mutant2sample)

        mutants2sample=next_round_mutants
        variant.append(sampled_mutant[0])
    print(nb_mutations)
    print(sampled_positions)
    variant= sorted(variant,key=lambda x:int(x[1:-1]))
    variant= ','.join(variant)
    assert variant not in training_variants, 'this would be highly unlikely but possible'
    if training_mutants:
        for mutant in variant.split(','):
            assert mutant not in training_mutants,'selected mutant should not be in training mutants for unobserved'
    Sequences.append(variant2seq(variant,av_gfp_WT))
    Variants.append(variant)
    Mutations.append(nb_mutations)
    return Sequences,Variants,Mutations

def random_sample():
    '''
    this is the function to randomluy
    :return:
    '''
    random.seed(42)
    train_idx = np.loadtxt(os.path.join('data','train_64_variants_avgfp.txt'),dtype=int)
    df= pd.read_csv(os.path.join('data','avgfp.csv'))

    train_df = df.iloc[train_idx]

    rejections = ",".join(train_df['variant'])
    training_mutants = rejections.split(',')
    print(f"training mutants: {len(training_mutants)}")
    training_mutants = list(np.unique(training_mutants))




    training_variants=list(train_df['variant'].apply(lambda x:sorted(x.split(','),
                                                                     key=lambda y:int(y[1:-1]))))
    all_mutants = get_all_singles(av_gfp_WT)


    unobserved_mutants =all_mutants
    for training_mutant in training_mutants:
        if training_mutant in unobserved_mutants:
            unobserved_mutants.remove(training_mutant)

    # Unobserved Mutants
    Variants =[]
    Mutations  =[]
    Label = []
    Sequences= [ ]
    for i in range(10):
        if i<5:
            Sequences,Variants,Mutation=add_variant(unobserved_mutants,5,Sequences,Variants,Mutations,training_variants,training_mutants)
        else:
            Sequences,Variants,Mutations=add_variant(unobserved_mutants, 10, Sequences, Variants, Mutations,training_variants,training_mutants)
        Label.append('Unobserved')

    for i in range(10):
        if i < 5:
            Sequences,Variants,Mutations=add_variant(training_mutants, 5, Sequences, Variants, Mutations,training_variants)
        else:
            Sequences,Variants,Mutations=add_variant(training_mutants, 10, Sequences, Variants, Mutations,training_variants)
        Label.append('Observed')

    df=pd.DataFrame(columns=['label','mutations','variant','sequence','seed'],data=np.array([Label,Mutations,Variants,Sequences,[42]*20]).T)
    df.to_csv('results/random_observed_unobserved/random_baseline.csv',index=False)

def verify_all_random_results():
    train_idx = np.loadtxt(os.path.join('data', 'train_64_variants_avgfp.txt'), dtype=int)
    df = pd.read_csv(os.path.join('data', 'avgfp.csv'))

    train_df = df.iloc[train_idx]
    rejections = ",".join(train_df['variant'])
    training_mutants = rejections.split(',')
    print(f"training mutants: {len(training_mutants)}")
    training_mutants = list(np.unique(training_mutants))

    tested_df= pd.read_csv('results/random_observed_unobserved/random_baseline.csv')

    for tested_variant,tested_label,tested_seq in zip(tested_df['variant'],tested_df['label'],tested_df['sequence']):
        Pos =  []
        new_seq =list(av_gfp_WT)
        for tested_mutant in tested_variant.split(','):
            wt_aa,pos,mut_aa= tested_mutant[0],int(tested_mutant[1:-1]),tested_mutant[-1]
            Pos.append(pos)

            if tested_label=='Unobserved':
                assert tested_mutant not in training_mutants
            if tested_label=='Observed':
                assert tested_mutant in  training_mutants

            assert list(av_gfp_WT)[pos]==wt_aa
            assert list(tested_seq)[pos]==mut_aa
            new_seq[pos]= mut_aa

        assert ''.join(new_seq)==tested_seq
        assert len(Pos)==len(list(np.unique(Pos)))
if __name__ == '__main__':
    verify_all_random_results()
    # random_sample()
    # check_constraints_all_mutants()
    # check_output_for_constraints()
    # preprocess_only_allow_training()
    # print_high_bin_variants()