import metl
import torch
import os
# the GFP wild-type sequence
wt ="SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQ" \
     "HDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKN" \
     "GIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
model, data_encoder = metl.get_from_checkpoint(os.path.join('src',"PEkeRuxb.pt"))
def seq2fitness_3d(mutants=None,WT=wt):
    '''
    3d sequence
    For avGFP, grab 1gfl_cm.pdb from https://github.com/samgelman/RosettaTL/tree/master/data/pdb_files
    :param mutants:
    :param WT:
    :return:
    '''
    # deal with edge case
    if type(mutants)==str:
        mutants=[mutants]
    if mutants==None:
        mutants = ["E3K,G102S"]
                    # "T36P,S203T,K207R",
                    # "V10A,D19G,F25S,E113V"]
    # some example GFP variants to compute the scores for

    encoded_variants = data_encoder.encode_variants(WT, mutants)

    # set model to eval mode
    model.eval()
    # no need to compute gradients for inference
    with torch.no_grad():
        predictions = model(torch.tensor(encoded_variants),pdb_fn=os.path.join('src','1gfl_cm.pdb'))
        

    # just return back a single scalar fitness prediction
    return float(predictions[0][0])

if __name__ == '__main__':
    predictions=seq2fitness_3d()
    print(predictions)