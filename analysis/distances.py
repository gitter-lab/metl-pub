
import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices
import time
import os
from src.utils import onehot,variant2sequence,av_gfp_WT
from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import cdist




# init blosum62 variants
blosum62 = substitution_matrices.load("BLOSUM62")
BLOSUM62_AA_libary=  list(blosum62.alphabet)
BLOSUMAA2Idx_MAPPING = {c: i for i, c in enumerate(BLOSUM62_AA_libary)}
BLOSUM62_max  = blosum62.max()

# create a pure numpy blosum62 matrix
# this is necessary for the vectorized function which runs 100x faster than direct string comparison
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

# helper function for vectorized blosum similarity
def onehot2index(onehot_encoding):
    idxs= []
    for row in onehot_encoding.reshape(-1,20):
        idxs.append(np.argmax(row))
    return np.array(idxs)

def blosum62_seq_similarity_strings(seq1,seq2):
    tot =0
    for aa1,aa2 in zip(seq1,seq2):
        tot+= blosum62[aa1,aa2]
    return tot

def calculate_normalized_dist_matrix(dist_matrix):
    WT_BLOSUM = blosum62_seq_similarity_strings(av_gfp_WT,av_gfp_WT)
    calculated_max = (BLOSUM62_max - 4) * 10 + WT_BLOSUM
    calculated_min = 0
    normalized_dist_matrix = 1 - ((dist_matrix - calculated_min) / (calculated_max - calculated_min))
    return normalized_dist_matrix

def calculated_dist_matrix_onehot(top_sequences_path,save_path):
    '''
    distance matrix onehot encoding
    :param top_sequences_path:the path to load in the top sequences from the simulated annealing runs
    :param save_path: where to save the distance matrix of the onehot encodings
    :return:
    '''
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['onehot'] = df['best_sequence'].apply(lambda x: onehot(x).flatten())
    X = np.vstack(df['onehot'].to_numpy())
    start= time.time()
    dist_matrix_onehot = cdist(X, X, metric='hamming')
    stop  = time.time()
    print(stop-start)
    np.save(os.path.join(save_path,'onehot_dist_matrix.npy'),dist_matrix_onehot)


def compute_distances(top_sequences_path,save_path):
    df = pd.read_csv(top_sequences_path)
    df['best_sequence'] = df['best_mutant'].apply(lambda x: variant2sequence(x))
    df['blosum_integer'] = df['best_sequence'].apply(lambda x: aa_seq_to_int_array(x))
    X = np.vstack(df['blosum_integer'].to_numpy())
    dist_matrix = pairwise_distances(X,metric=lambda u,v: vectorized_blosum_similarity(u,v,blosum_matrix=blosum62_arr))
    np.save(os.path.join(save_path,'vectorized_blosum_similarity_distance_matrix.npy'),dist_matrix)
