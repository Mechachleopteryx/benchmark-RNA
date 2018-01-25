#!/usr/bin/env python3
import sys
import os
import shlex
import subprocess
from subprocess import Popen, PIPE, STDOUT
import re
import copy
import argparse
import glob


#warnings
def printWarningMsg(msg):
	print("[Warning] " + msg)


# check if file exists and is not empty
def checkIfFile(pathToFile):
	if not(os.path.exists(pathToFile) and os.path.getsize(pathToFile) > 0):
		return False
	return True

# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	#~ p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
	p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
	return p

# find files with a regex
def getFiles(pathToFiles, name): #for instance name can be "*.txt"
	os.chdir(pathToFiles)
	listFiles = []
	for files in glob.glob(name):
		listFiles.append(files)
	return listFiles

# read simulation
def simulateReads(covSR, covLR, skipped, abund):
	for sizeSkipped in skipped:
		for relAbund in  abund:
			suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
			# simulation
			cmdSimul = "./ES_simulation " + str(sizeSkipped) + " " + str(relAbund) + " " + str(suffix) + " " +  str(covSR) + " " + str(covLR)
			cmdSimul = subprocessLauncher(cmdSimul)

# return number of reads in a fasta
def getFileReadNumber(fileName):
	cmdGrep = """grep ">" -c """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return int(val.decode('ascii'))

# return the sequence of an errorless isoform
def getPerfectSequence(fileName):
	cmdGrep = """grep "[ACGT]" -m 1 """ + fileName
	seq = subprocess.check_output(['bash','-c', cmdGrep])
	return seq.decode('ascii')


# associate to isoform type the headers of the reference file
def makeReferenceHeadersList(currentDirectory, skipped, abund):
	listFilesPerfect = getFiles(currentDirectory, "perfect*_size_" + skipped + "_abund_" + abund + ".fa")
	refIsoformTypesToCounts = dict()
	refIsoformTypesToSeq = dict()
	for fileP in listFilesPerfect:
		typeIsoform = fileP.split("_")[2]
		readNb = getFileReadNumber(currentDirectory + "/" + fileP) #get nb of reads
		#~ noIsoform = range(readNb)
		headers = [typeIsoform + str(x) for x in range(readNb)]
		refIsoformTypesToCounts[typeIsoform] = headers
		perfectSeq = getPerfectSequence(currentDirectory + "/" + fileP)
		refIsoformTypesToSeq[typeIsoform] = perfectSeq.rstrip()
	return (listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq)

# align triplets of reads
def msa(suffix, msaType):
	print("****************************", msaType)
	if msaType == "msa_isoform":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa"
	elif msaType == "msa_exon":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa -c exon"
	elif msaType == "msa_sparc":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa -s True"
	p = subprocessLauncher(cmdMSA)

# headers of corrected reads file
def getCorrectedHeaders(fileName):
	cmdGrep = """grep ">" """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return val.decode('ascii').split("\n")[:-1] #last item is empty


# get consensus sequence from corrected reads file
def getCorrectedSequence(fileName):
	cmdGrep = """grep ">" -v -m 1 """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return val.decode('ascii').rstrip().upper() #last item is empty


# compare headers coming from an isoform in ref file to headers attributed to an isoform in correction file
def compareRefAndCorrectedHeaders(refIsoformTypesToCounts, correcIsoformTypesToCounts):
	count = dict()
	for typeIsoform in refIsoformTypesToCounts.keys():
		count[typeIsoform] = dict()
		for correcIsoform in correcIsoformTypesToCounts.keys():
			for headersRef in refIsoformTypesToCounts[typeIsoform]:
				for headersCor in correcIsoformTypesToCounts[correcIsoform]:
					if headersRef == headersCor:
						if correcIsoform in count[typeIsoform].keys():
							count[typeIsoform][correcIsoform] += 1
						else:
							count[typeIsoform][correcIsoform] = 1
	return count #for instance {'exclusion': {'exclusion': 5}, 'inclusion': {'inclusion': 5}}
	

# associate to isoform type the headers of the corrected file
def makeCorrectedHeadersList(resultDirectory, currentDirectory, skipped, abund, suffix, refIsoformTypesToCounts):
	listFilesCorrected = getFiles(resultDirectory, "corrected_by_MSA*.fa")
	correcIsoformTypesToCounts = dict()
	correcIsoformTypesToSeq = dict()
	for fileC in listFilesCorrected:
		correctedSequence = getCorrectedSequence(fileC)
		listHeaders = getCorrectedHeaders(fileC)
		listHeaders = [x.split("_")[0][1:] + x.split("_")[1].split(' ')[0] for x in listHeaders] #For instance ['exclusion0', 'exclusion1', 'exclusion2', 'exclusion3', 'exclusion4']
		listIsoforms = [x.split("_")[0][:-1] for x in listHeaders] 
		for isoformType in set(listIsoforms): #unique types of isoforms
			correcIsoformTypesToCounts[isoformType] = listHeaders
			correcIsoformTypesToSeq[isoformType] = correctedSequence
	return(correcIsoformTypesToCounts, correcIsoformTypesToSeq)

#compute ratio of isoforms representations in ref and corrected files
def computeRatioIsoforms(refIsoformTypesToCounts, correcIsoformTypesToCounts, currentDirectory, suffix):
	counts = compareRefAndCorrectedHeaders(refIsoformTypesToCounts, correcIsoformTypesToCounts)
	confusionName = currentDirectory + "/matrix_confusion" + suffix + ".txt"
	outConf = open(confusionName, 'w')
	outConf.write("reference correction ratio\n")
	for ref in counts:
		for ref2 in counts:
			if ref2 in counts[ref].keys():
				ratio = counts[ref][ref2] * 1.0 / len(refIsoformTypesToCounts[ref]) if len(refIsoformTypesToCounts[ref]) != 0 else 0
				outConf.write(ref + " " + ref2 + " " + str(ratio) + "\n")
			else:
				outConf.write(ref + " " + ref2 + " 0\n")
	outConf.close()
	return confusionName




#R functions
def printConfusionMatrix(currentDirectory, corrector, confusionFile, suffix):
	Rcmd = "Rscript " + currentDirectory + "/matrice_confusion.R " + confusionFile + " " + corrector  + suffix + " " + currentDirectory
	subprocessLauncher(Rcmd)


def computeResults(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts, outDir="/home/marchet/detection-consensus-isoform/results"):
	msa(suffix, correc)
	correcIsoformTypesToCounts, correcIsoformTypesToSeq = makeCorrectedHeadersList(outDir, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts)
	confusionFile = computeRatioIsoforms(refIsoformTypesToCounts, correcIsoformTypesToCounts, currentDirectory, suffix)
	printConfusionMatrix(currentDirectory, correc, confusionFile, suffix)


def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))
	covSR = 1
	covLR = 20
	#~ skipped = [50,100]
	#~ abund = [50,75,90]
	abund = [10]
	skipped = [100]
	#~ abund = [50]
	skippedS = [str(r) for r in skipped]
	abundS = [str(r) for r in abund]

	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	parser = argparse.ArgumentParser()
	parser.add_argument('-output', nargs='?', type=str, action="store", dest="outputDirPath", help="Name for output directory", default=None)
	#~ parser.add_argument('-corrector', type=str, action="store", dest="correctors", help="A particular corrector to be used", default=None)
	parser.add_argument('-coverageLR', nargs='?', type=int, action="store", dest="covLR", help="Coverage for LR simulation (default 10)", default=10)
	# get options for this run
	args = parser.parse_args()
	outputDirPath = args.outputDirPath
	covLR = args.covLR

	correctors = ["msa_isoform", "msa_exon"]
	#~ correctors = ["msa_isoform"]
	#~ correctors = ["msa_exon"]
	if not outputDirPath is None:
		if not os.path.exists(outputDirPath):
			os.mkdir(outputDirPath)
		else:
			printWarningMsg(outputDirPath+ " directory already exists, we will use it.")
			try:
				cmdRm = "(cd " + outputDirPath + " && rm *)"
				subprocess.check_output(['bash','-c', cmdRm])
			except subprocess.CalledProcessError:
				pass
	simulateReads(covSR, covLR, skipped, abund)
	for correc in correctors:
		for skippedExon in skippedS:
			for abundanceMajor in abundS:
				listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq = makeReferenceHeadersList(currentDirectory, str(skippedExon), str(abundanceMajor))
				suffix = "_size_" + str(skippedExon) + "_abund_" + str(abundanceMajor)
			
				computeResults(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts)

if __name__ == '__main__':
	main()
