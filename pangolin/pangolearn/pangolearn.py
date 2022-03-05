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
import joblib
from Bio import SeqIO
import sys
import multiprocessing as mp
import os
import warnings

warnings.filterwarnings("ignore")

def findReferenceSeq(referenceFile):
	currentSeq = ""

	with open(referenceFile) as f:
		for line in f:
			if ">" not in line:
				currentSeq = currentSeq + line.strip()

	f.close()
	return currentSeq

# function for handling weird sequence characters
def clean(referenceSeq,x, loc):
	x = x.upper()

	if x == 'T' or x == 'A' or x == 'G' or x == 'C' or x == '-':
		return x, False

	if x == 'U':
		return 'T', False

	# replace ambiguity with the reference seq value

	if referenceSeq[loc] == x:
		return referenceSeq[loc], False

	return referenceSeq[loc], True

# generates data line
def encodeSeq(referenceSeq,seq, indiciesToKeep):
	dataLine = []
	imputed = 0
	nonimputed = 0

	for i in indiciesToKeep:
		if i < len(seq):
			cleaned, isImputed = clean(referenceSeq,seq[i], i)

			dataLine.extend(cleaned)

			if(isImputed):
				imputed = imputed + 1
			else:
				nonimputed = nonimputed + 1

	score = 0

	if nonimputed != 0 and imputed < nonimputed:
		score = 1 - (imputed/nonimputed)

	return dataLine, score

# reads in the two data files
def readInAndFormatData(referenceSeq, imputationScores,sequencesFile, indiciesToKeep, blockSize=1000):
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

						finalSeq, imputationScore = encodeSeq(referenceSeq,currentSeq, indiciesToKeep)
						seqList.append(finalSeq)
						currentSeq = ""

						imputationScores[seqid] = imputationScore

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
		finalSeq, imputationScore = encodeSeq(referenceSeq,currentSeq, indiciesToKeep)
		seqList.append(finalSeq)
		imputationScores[seqid] = imputationScore

	yield idList, seqList


def file_writer(outfile, q):
	"""listens for messages on the q, writes to files. """
	f = open(outfile, "w")
	f.write("taxon,prediction,score,imputation_score,non_zero_ids,non_zero_scores,designated\n")
	while 1:
		message = q.get()
		if message == 'kill':
			break
		f.write(message)
	f.close()


def assign_to_chunk(idList, seqList, refRow, referenceId, loaded_model, imputationScores, columns, categories, q):
	print("Processing block of {} sequences {}".format(
		len(seqList), datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
	))

	rows = [[r == c for r in row for c in categories] for row in seqList]
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

		nonZeroIds = []
		nonZeroScores = []

		# get the max probability score and its assosciated index
		for i in range(len(predictions[index])):
			if predictions[index][i] > maxScore:
				maxScore = predictions[index][i]
				maxIndex = i

				nonZeroScores.append(predictions[index][i])
				nonZeroIds.append(loaded_model.classes_[i])

		score = maxScore
		prediction = loaded_model.classes_[maxIndex]
		seqId = idList[index]

		nonZeroIds = ";".join(nonZeroIds)
		nonZeroScores = ';'.join(str(x) for x in nonZeroScores)

		if seqId != referenceId:
			res = seqId + "," + prediction + "," + str(score) + "," + str(
				imputationScores[seqId]) + "," + nonZeroIds + "," + nonZeroScores + "," + "\n"
			q.put(res)
	return 0


def assign_lineage(header_file,model_file,reference_file,sequences_file,outfile,threads=1):
	print("Running pangoLEARN assignment")
	dirname = os.path.dirname(__file__)

	referenceSeq = ""
	referenceId = "reference"

	imputationScores = dict()

	# setup manager
	manager = mp.Manager()
	q = manager.Queue()
	pool = mp.Pool(threads)

	if threads > 1:
		pool.apply_async(file_writer, (outfile, q))

	records = 0
	for record in SeqIO.parse(sequences_file,"fasta"):
		records +=1

	# fire off workers
	jobs = []

	# run batches
	if records ==0:
		q.put("kill")
		print("No sequences to assign with pangoLEARN.")
	else:

		# loading the list of headers the model needs.
		model_headers = joblib.load(header_file)
		indiciesToKeep = model_headers[1:]

		referenceSeq = findReferenceSeq(reference_file)
		# possible nucleotide symbols
		categories = ['-','A', 'C', 'G', 'T']
		columns = [f"{i}_{c}" for i in indiciesToKeep for c in categories]



		rs, score = encodeSeq(referenceSeq, referenceSeq, indiciesToKeep)

		refRow = [r==c for r in rs for c in categories]

		print("Loading model " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
		loaded_model = joblib.load(model_file)
		print("Finished loading model " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))

		# write predictions to a file
		for idList, seqList in readInAndFormatData(referenceSeq,imputationScores,sequences_file, indiciesToKeep):
			job = pool.apply_async(assign_to_chunk, (idList, seqList, refRow, referenceId, loaded_model,
													 imputationScores, columns, categories, q))
			jobs.append(job)

	# collect results from the workers through the pool result queue
	for job in jobs:
		job.get()

	# now we are done, kill the listener
	q.put('kill')

	if threads == 1:
		pool.apply_async(file_writer, (outfile, q))

	pool.close()
	pool.join()

	print("Complete " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))