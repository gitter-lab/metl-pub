import metl
import torch


# the GFP wild-type sequence
wt = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQ" \
     "HDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKN" \
     "GIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
model, data_encoder = metl.get_from_checkpoint("YoQkzoLD.pt")
def seq2fitness(mutants=None,WT=wt):
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
        predictions = model(torch.tensor(encoded_variants))

    return float(predictions[0][0])

if __name__ == '__main__':
    predictions=seq2fitness()
    print(predictions)
