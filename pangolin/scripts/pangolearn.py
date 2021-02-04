#!/usr/bin/env python3

import pandas as pd
import numpy as np
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
import os

def parse_args():
	parser = argparse.ArgumentParser(description='pangoLEARN.')
	parser.add_argument("--header-file", action="store", type=str, dest="header_file")
	parser.add_argument("--model-file", action="store", type=str, dest="model_file")
	parser.add_argument("--fasta", action="store", type=str, dest="sequences_file")
	parser.add_argument("--reference-file",action="store",type=str,dest="reference_file")
	parser.add_argument("-o","--outfile", action="store", type=str, dest="outfile")
	return parser.parse_args()


args = parse_args()

dirname = os.path.dirname(__file__)
referenceFile = args.reference_file

referenceSeq = ""
referenceId = "reference"

def findReferenceSeq():
	currentSeq = ""

	with open(referenceFile) as f:
		for line in f:
			if ">" not in line:
				currentSeq = currentSeq + line.strip()

	f.close()
	return currentSeq


# function for handling weird sequence characters
def clean(x, loc):
	x = x.upper()

	if x == 'T' or x == 'A' or x == 'G' or x == 'C' or x == '-':
		return x

	if x == 'U':
		return 'T'

	# replace ambiguity with the reference seq value
	return referenceSeq[loc]


# generates data line
def encodeSeq(seq, indiciesToKeep):
	dataLine = []
	for i in indiciesToKeep:
		if i < len(seq):
			dataLine.extend(clean(seq[i], i))

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
						currentSeq = ""

					# this is a fasta line designating an id, but we don't want to keep the >
					seqid = line.strip('>')

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
model_headers = joblib.load(args.header_file)
indiciesToKeep = model_headers[1:]

referenceSeq = findReferenceSeq()
# possible nucleotide symbols
categories = ['-','A', 'C', 'G', 'T']
columns = [f"{i}_{c}" for i in indiciesToKeep for c in categories]
refRow = [r==c for r in encodeSeq(referenceSeq, indiciesToKeep) for c in categories]

print("loading model " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
loaded_model = joblib.load(args.model_file)

# write predictions to a file
f = open(args.outfile, "w")

for idList, seqList in readInAndFormatData(args.sequences_file, indiciesToKeep):
	print("processing block of {} sequences {}".format(
		len(seqList), datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
	))

	rows = [[r==c for r in row for c in categories] for row in seqList]
	# the reference seq must be added to everry block to make sure that the 
	# spots in the reference have Ns are in the dataframe to guarentee that 
	# the correct number of columns is created when get_dummies is called
	rows.append(refRow)
	idList.append(referenceId)

	# create a data from from the seqList
	d = np.array(rows, np.uint8)
	df = pd.DataFrame(d, columns=columns)

	predictions = loaded_model.predict_proba(df)

	for index in range(len(predictions)):

		maxScore = 0
		maxIndex = -1

		# get the max probability score and its assosciated index
		for i in range(len(predictions[index])):
			if predictions[index][i] > maxScore:
				maxScore = predictions[index][i]
				maxIndex = i

		score = maxScore
		prediction = loaded_model.classes_[maxIndex]
		seqId = idList[index]

		if seqId != referenceId:
			f.write(seqId + "," + prediction + "," + str(score) + "\n")

f.close()

print("complete " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
