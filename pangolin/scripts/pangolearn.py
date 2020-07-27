import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from datetime import datetime
import argparse
import joblib
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='pangoLEARN.')
    parser.add_argument("--header-file", action="store", type=str, dest="header_file")
    parser.add_argument("--model-file", action="store", type=str, dest="model_file")
    parser.add_argument("--fasta", action="store", type=str, dest="sequences_file")
    parser.add_argument("-o","--outfile", action="store", type=str, dest="outfile")
    return parser.parse_args()


args = parse_args()


# small class to store vector objects
class VectorObject:
    def __init__(self, vector):
        self.vector = vector

    def equals(self, vector):
        for i in range(len(vector.vector)):
            if vector.vector[i] != self.vector[i]:
                return False

        return True


# produces the proper one-hot encoding for a particular genomic site
def getOneHotEncoding(char):
    alphabet = 'ATCGN-'
    v = [0] * len(alphabet)
    if char in alphabet:
        v[alphabet.index(char)] = 1
        return VectorObject(v)
    else:
        # default to N
        return VectorObject([0, 0, 0, 0, 1, 0])


# converts sequence to flattened one-hot vector
def encodeSeq(seq, indiciesToKeep):
    dataLine = []
    for i in indiciesToKeep:
        if i < len(seq):
            dataLine.extend(getOneHotEncoding(seq[i]).vector)
    return dataLine


# reads in the two data files
def readInAndFormatData(sequencesFile, indiciesToKeep, blockSize=1000):
    idList = []
    seqList = []

    # open sequencesFile, which is the first snakemake input file
    with open(sequencesFile) as f:
        currentSeq = ""

        # go line by line through the file, collecting a list of the sequences
        for line in f:
            # if the line isn't the header line
            if "taxon,lineage" not in line:
                line = line.strip()

                if ">" in line:
                    # starting new entry, gotta save the old one
                    if currentSeq:
                        # yield sequence as one-hot encoded vector
                        idList.append(seqid)
                        seqList.append(encodeSeq(currentSeq, indiciesToKeep))

                    # this is a fasta line designating an id, but we don't want to keep the >
                    seqid = line.strip('>')
                    currentSeq = ""

                else:
                    currentSeq = currentSeq + line

            if len(seqList) == blockSize:
                yield idList, seqList
                idList = []
                seqList = []

        # gotta get the last one
        idList.append(seqid)
        seqList.append(encodeSeq(currentSeq, indiciesToKeep))

    yield idList, seqList



# loading the list of headers the model needs.
# example: '29657-A', '29657-T', '29657-C', '29657-G', '29657-N', '29657-gap'
model_headers = joblib.load(args.header_file)
indiciesToKeep = []

# by cycling through model_headers, get which column indicies we need to keep in the test data
for h in model_headers:
    if h != "lineage" and "-A" in h:
        # -1 because in the training data the 0 position is the lineage
        index = int(h.split("-")[0]) - 1
        indiciesToKeep.append(index)


print("loading model " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
loaded_model = joblib.load(args.model_file)

print("generating predictions " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
# write predictions to a file
f = open(args.outfile, "w")

for idList, seqList in readInAndFormatData(args.sequences_file, indiciesToKeep):
    print("processing block of {} sequences {}".format(
        len(seqList), datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
    ))
    df = pd.DataFrame(seqList, columns=model_headers[1:])
    predictions = loaded_model.predict_proba(df)

    for row, pred in enumerate(predictions):
        # find index of maximum prediction score
        intermed = [(pr, idx) for idx, pr in enumerate(pred)]
        intermed.sort(reverse=True)
        maxScore, maxIndex = intermed[0]

        prediction = loaded_model.classes_[maxIndex]
        f.write("{},{},{}\n".format(idList[row], prediction, maxScore))

f.close()

print("complete " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
