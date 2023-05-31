import numpy as np
AAs = 'ACDEFGHIKLMNPQRSTVWY'

CHARS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
C2I_MAPPING = {c: i for i, c in enumerate(CHARS)}
av_gfp_WT= 'SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLT' \
    'LKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMP' \
    'EGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGN' \
    'ILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQ' \
    'QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGI' \
    'THGMDELYK'
def variant2sequence(variants,WT=av_gfp_WT):
    WT_list= list(WT)
    for variant in variants.split(','):
        assert WT_list[int(variant[1:-1])] == variant[0], 'wt and variant wt must be equal'
        WT_list[int(variant[1:-1])] = variant[-1]
    return "".join(WT_list)

def onehot(sequence):
    '''
    take a sequence and return onehot encoding
    :param sequence:
    :return:
    '''
    X= np.zeros((len(sequence),len(CHARS)))
    for i,s in enumerate(sequence):
        X[i][C2I_MAPPING[s]]= 1
    return X

def onehot2sequence(onehot):
    reshaped = onehot.reshape(-1, 20)
    sequence  = ''
    for row in reshaped:
        sequence+= CHARS[np.argmax(row)]
    return sequence



