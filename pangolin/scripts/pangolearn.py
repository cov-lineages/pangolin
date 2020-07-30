#!/usr/bin/env python3

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

dataList = []
tempDataLines = []
idList = []

def parse_args():
    parser = argparse.ArgumentParser(description='pangoLEARN.')
    parser.add_argument("--header-file", action="store", type=str, dest="header_file")
    parser.add_argument("--model-file", action="store", type=str, dest="model_file")
    parser.add_argument("--fasta", action="store", type=str, dest="sequences_file")
    parser.add_argument("-o","--outfile", action="store", type=str, dest="outfile")
    return parser.parse_args()


args = parse_args()

sequencesFile = args.sequences_file
modelFile = args.model_file
headerFile = args.header_file

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
	# ATCGN-
	if char == "A":
		return VectorObject([1, 0, 0, 0, 0, 0])
	elif char == "T":
		return VectorObject([0, 1, 0, 0, 0, 0])
	elif char == "C":
		return VectorObject([0, 0, 1, 0, 0, 0])
	elif char == "G":
		return VectorObject([0, 0, 0, 1, 0, 0])
	elif char == "-":
		return VectorObject([0, 0, 0, 0, 0, 1])
	else:
		return VectorObject([0, 0, 0, 0, 1, 0])


# reads in the two data files
def readInAndFormatData():
	# open sequencesFile, which is the first snakemake input file
	with open(sequencesFile) as f:
		currentSeq = ""

		# go line by line through the file, collecting a list of the sequences
		for line in f:
			# if the line isn't the header line
			if "taxon,lineage" not in line:
				line = line.strip()

				if ">" in line:
					# this is a fasta line designating an id, but we don't want to keep the >
					idList.append(line[1:])

					# starting new entry, gotta save the old one
					if currentSeq:
						tempDataLines.append(currentSeq)

					currentSeq = ""

				else:
					currentSeq = currentSeq + line

		# gotta get the last one
		if  currentSeq:
			tempDataLines.append(currentSeq)

	# get a list of viral genomes, represented as a list of one-hot encoded vectors
	for line in tempDataLines:
		dataLine = []
		for index in range(len(line)):
				char = line[index]

				dataLine.append(getOneHotEncoding(char))
				
		dataList.append(dataLine)

	# close the file
	f.close()


# remove the indicies representing genomic locations which aren't used by the model
def removeIndices(headersFile):
	# will hold the final data from which the model will make predictions
	finalList = []
	# list of headers representing the genomic locations of each entry
	headers = list(range(len(dataList[0])))
	# list which will hold the indicies of interest
	indiciesToKeep = []

	# loading the list of headers the model needs.
	# example: '29657-A', '29657-T', '29657-C', '29657-G', '29657-N', '29657-gap'
	model_headers = joblib.load(headersFile)

	# by cycling through model_headers, get which column indicies we need to keep in the test data
	for h in model_headers:
		if h != "lineage" and "-A" in h:
			# -1 because in the training data the 0 position is the lineage
			index = int(h.split("-")[0]) - 1
			indiciesToKeep.append(index)

	# for each entry in dataList, remove the irrelevant columns
	while len(dataList) > 0:
		line = dataList.pop(0)

		finalLine = []

		for index in range(len(line)):
			if index in indiciesToKeep:
				finalLine.extend(line[index].vector)

		finalList.append(finalLine)

	# return the final data, along with the relevant headers 
	# model_headers[1:] because the first is just "lineage" which we won't have for this test data
	return finalList, model_headers[1:]


print("reading in data " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"));

readInAndFormatData()

print("removing unnecessary columns " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"));

dataList, headers = removeIndices(headerFile)

df = pd.DataFrame(dataList, columns=headers)

print("loading model " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"));

loaded_model = joblib.load(modelFile)

print("generating predictions " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"));

predictions = loaded_model.predict_proba(df)

# write predictions to a file
f = open(args.outfile, "w")
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

	f.write(seqId + "," + prediction + "," + str(score) + "\n")

f.close()
