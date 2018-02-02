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

# to return if an error makes the run impossible
def dieToFatalError (msg):
  print("[FATAL ERROR] " + msg)
  sys.exit(1);


# check if file exists and is not empty
def checkIfFile(pathToFile):
	if not(os.path.exists(pathToFile) and os.path.getsize(pathToFile) > 0):
		return False
	return True

# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
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
def simulateReads(skipped, abund, coverage, EStype, currentDirectory, errorRate):
	covSR = 10
	for error in errorRate:
		for covLR in coverage:
			for sizeSkipped in skipped:
				for relAbund in abund:
					suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund) + "_cov_" + str(covLR) + "_err_" + str(error)
					# simulation
					if checkIfFile(currentDirectory + "/simulatedLR"+ suffix +".fa" ):
						cmdRm = "rm " + currentDirectory + "/simulatedLR"+ suffix +".fa"
						subprocess.check_output(['bash','-c', cmdRm])
					if EStype == "ES":
						cmdSimul = currentDirectory + "/ES_simulation " + str(sizeSkipped) + " " + str(relAbund) + " " + str(suffix) + " " +  str(covSR) + " " + str(covLR) + " " + str(error)
					elif EStype == "MES": #todo make error rate a parameter for this one
						cmdSimul =  currentDirectory + "/MES_simulation " + str(sizeSkipped) + " " + str(relAbund) + " " + str(suffix) +  " " + str(covLR)
					elif EStype == "alt":
						cmdSimul =  currentDirectory + "/AltSE_simulation " + str(sizeSkipped) + " " + str(relAbund) + " " + str(suffix) + " " +  str(covSR) + " " + str(covLR) + " " + str(error)
					cmdSimul = subprocessLauncher(cmdSimul)
					checkReadFiles(currentDirectory + "/simulatedLR"+ suffix +".fa")


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

	
def getPerfectSequenceLength(fileName):
	cmdWc = """grep "[ACGT]" -m 1 """ + fileName + "| wc"
	seq = subprocess.check_output(['bash','-c', cmdWc])
	return int(seq.decode('ascii').split(" ")[-1].rstrip())


# associate to isoform type the headers of the reference file
def makeReferenceHeadersList(currentDirectory, skipped, abund, cov, error):
	listFilesPerfect = getFiles(currentDirectory, "perfect*_size_" + skipped + "_abund_" + abund + "_cov_" + cov  + "_err_" + str(error) + ".fa")
	print(listFilesPerfect)
	refIsoformTypesToCounts = dict()
	refIsoformTypesToSeq = dict()
	for fileP in listFilesPerfect:
		typeIsoform = fileP.split("_")[2]
		readNb = getFileReadNumber(currentDirectory + "/" + fileP) #get nb of reads
		#~ noIsoform = range(readNb)
		headers = [typeIsoform + str(x) for x in range(readNb)]
		print(headers)
		refIsoformTypesToCounts[typeIsoform] = headers
		perfectSeq = getPerfectSequence(currentDirectory + "/" + fileP)
		refIsoformTypesToSeq[typeIsoform] = perfectSeq.rstrip()
	return (listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq)

# align triplets of reads
def msa(suffix, msaType,outDir = "/home/marchet/detection-consensus-isoform/results"):
	print("Launch", msaType)
	if msaType == "msa_isoform":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa -c isoform"
	elif msaType == "msa_exon":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa "
	elif msaType == "msa_sparc":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa -s True"
	elif msaType == "msa_both":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py -r simulatedLR" + suffix + ".fa -c both"
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
# todo files with reference correction ratio soft coverage errorrate 
def computeRatioIsoforms(refIsoformTypesToCounts, correcIsoformTypesToCounts, currentDirectory, suffix, soft, coverage, errorRate):
	counts = compareRefAndCorrectedHeaders(refIsoformTypesToCounts, correcIsoformTypesToCounts)
	confusionName = currentDirectory + "/matrix_confusion" + suffix +  ".txt"
	outConf = open(confusionName, 'w')
	outConf.write("reference correction ratio\n")
	isCorrect = True
	ratios = dict()
	for ref in counts:
		for ref2 in counts:
			if ref2 in counts[ref].keys():
				ratio = counts[ref][ref2] * 1.0 / len(refIsoformTypesToCounts[ref]) if len(refIsoformTypesToCounts[ref]) != 0 else 0
				if ratio != 1:
					isCorrect = False
				if ref in ratios.keys():
					if ref2 in ratios[ref].keys():
						ratios[ref][ref2].append(ratio)
					else:
						ratios[ref][ref2] = [ratio]
				else:
					ratios[ref] = dict()
					ratios[ref][ref2] = [ratio]
	for ref in ratios.keys():
		for ref2 in ratios.keys():
			if ref2 in ratios[ref].keys():
				meanRatio = sum(ratios[ref][ref2])/len(ratios[ref][ref2]) if len(ratios[ref][ref2]) > 0 else 0
				outConf.write(ref + " " + ref2 + " " + str(meanRatio) + "\n")
				cmdEcho = "echo " + ref + " " + ref2 + " " + str(meanRatio) + " " + soft + " " + str(coverage) + " " + str(errorRate) + " >> " + currentDirectory + "/all_confusion_matrix.txt"
				subprocess.check_output(['bash','-c', cmdEcho])
			else:
				outConf.write(ref + " " + ref2 + " 0\n")
				cmdEcho = "echo " + ref + " " + ref2 + " 0 " + soft + " " + str(coverage) + " " + str(errorRate) + " >> " + currentDirectory + "/all_confusion_matrix.txt"
				subprocess.check_output(['bash','-c', cmdEcho])

	outConf.close()
	return confusionName, isCorrect




def alignOnRefMsa(soft, skipped, abund, currentDirectory, resultDirectory, cov, error):
	suffix = "_size_" + str(skipped) + "_abund_" + str(abund) + "_cov_" + str(cov)  + "_err_" + str(error)
	listFileNames = getFiles(resultDirectory, "corrected_by_MSA*.fa")
	for fileC in listFileNames:
		isoform = getCorrectedHeaders(resultDirectory + "/" + fileC)[0].split("_")[0][1:]
		cmdGrep = "grep "+ isoform + " " + currentDirectory + "/refSequences" + suffix + ".fa -A 1 > " + currentDirectory + "/refSequences" + isoform + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdGrep])
		cmdGrep = "grep "+ isoform + " " + fileC + " -A 1 > toalign.fa"
		subprocess.check_output(['bash','-c', cmdGrep])
		samFile = open(currentDirectory + "/results" + isoform + soft + suffix + ".sam", 'w')
		cmdAlign = "/home/marchet/bin/Complete-Striped-Smith-Waterman-Library/src/ssw_test " + currentDirectory +"/refSequences" + isoform + suffix + ".fa toalign.fa -c -s"
		p = subprocessLauncher(cmdAlign, samFile)
		samFile.close()
		cmdCp = "cp " + resultDirectory + "/" +  fileC + " " + currentDirectory + "/corrected_reads_by_" + soft + "_" + isoform + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdCp])




def readSam(soft, suffix, isoformType, currentDirectory):
	blockResults = dict()
	lenResults = dict()
	pathSam = currentDirectory + "/results" +isoformType + soft + suffix + ".sam"
	if os.path.exists(pathSam) and os.path.getsize(pathSam) > 0:
		samFile = open(pathSam, 'r')
		readsSize = []
		lines = samFile.readlines()
		queries = dict()
		for line in lines:
			line = line.rstrip().split('\t')
			query = line[0]
			target = line[2]
			cigar = line[5]
			start = int(line[3]) - 1
			length = len(line[9])
			seq = line[9]
			readsSize.append(length)
			blocks = re.compile("[0-9]+").split(cigar)[1:]
			resultAln = re.compile("[A-Z]|=").split(cigar)[:-1]
			alnLength = 0
			gapsLength = 0
			queries[query] = seq
			if len(blocks) == 1 and len(resultAln) == 1 and blocks[0] == '=': #aligned in one block
				blockResults[query] = {1:target}
				alnLength = int(resultAln[0])
			else:
				if query not in blockResults:
					blockResults[query] = {len(blocks):target}
				else:
					if 1 not in blockResults[query]:
						blockResults[query][len(blocks)] = target
				
				for b,r in zip(blocks, resultAln):
					
					if b != 'D' and b!= 'I':
						alnLength += int(r)
					else:
						gapsLength += int(r)
			if query not in lenResults:
				lenResults[query] = {target:[length, alnLength -start, gapsLength]}
			else:
				lenResults[query][target] = [length, alnLength -start, gapsLength]
		return ( start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries)
	else:
		return (None,) * 8



def computeResultsRecallPrecision(corrector, skipped, abund, currentDirectory, soft, refIsoformTypesToSeq, outSize, cov, error):
	suffix = "_size_" + str(skipped) + "_abund_" + str(abund) + "_cov_" + str(cov)  + "_err_" + str(error)
	cmdFile = "> precision_tmp.txt"
	subprocess.check_output(['bash','-c', cmdFile])
	cmdFile = "> recall_tmp.txt"
	subprocess.check_output(['bash','-c', cmdFile])
	cmdFile = "> correct_base_rate_tmp.txt"
	subprocess.check_output(['bash','-c', cmdFile])
	#~ expectedLengths = getExpectedLength(currentDirectory, suffix, isofType)
	#~ ratioLen = []
	for isofType in refIsoformTypesToSeq:
		expectedLengths = getPerfectSequenceLength(currentDirectory + "/perfect_reads_" + isofType + suffix + ".fa")
							

		start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries = readSam(soft, suffix, isofType, currentDirectory)
		#~ meanSizes[isofType] = {"realSize" : [], "alignedSize" : []}
		#~ for querySeq, aln in blockResults.items():
			#~ meanSizes[isofType]["realSize"].append(lenResults[querySeq][isofType][0])
			#~ meanSizes[isofType ]["alignedSize"].append(lenResults[querySeq][isofType][1])
		meanReadsSize = round(sum(readsSize)*1.0/len(readsSize),2) if len(readsSize) > 0 else 0
		ratioLen = round(meanReadsSize*100/expectedLengths,2)
		outSize.write(soft + " " + str(ratioLen) + " " +  str(skipped) + " "+ str(abund) +"\n")

		#~ ratioLenE = round(meanExclusionCorrectedSize*100/expectedLengths["exclusion"],2)
		
		cmdHead = "head -2 " + currentDirectory + "/corrected_reads_by_" + soft + "_" + isofType + suffix + ".fa > " + currentDirectory + "/corrected.fa"
		subprocess.check_output(['bash','-c', cmdHead])
		cmdHead = "head -2 " + currentDirectory + "/uncorrected_reads_"  + isofType + suffix + ".fa > " + currentDirectory + "/uncorrected.fa"
		subprocess.check_output(['bash','-c', cmdHead])
		cmdHead = "head -2 " + currentDirectory + "/perfect_reads_" + isofType + suffix + ".fa > " + currentDirectory + "/perfect.fa"
		subprocess.check_output(['bash','-c', cmdHead])
		cmdBench = "python3 " + currentDirectory + "/benchmark-long-read-correction/benchmark.py -c " + currentDirectory + "/corrected.fa -u " + currentDirectory + "/uncorrected.fa -r " + currentDirectory + "/perfect.fa -o " + currentDirectory
		subprocessLauncher(cmdBench)

		cmdSed = "sed -i 's/unknown/" + soft + "/g' " + currentDirectory + "/precision.txt"
		subprocess.check_output(['bash','-c', cmdSed])
		cmdCat = "grep " + soft + " " + currentDirectory + "/precision.txt >> " + currentDirectory + "/precision_tmp.txt"
		subprocess.check_output(['bash','-c', cmdCat])
		
		cmdSed = "sed -i 's/unknown/" + soft + "/g' " + currentDirectory + "/recall.txt"
		subprocess.check_output(['bash','-c', cmdSed])
		cmdCat = "grep " + soft +" "  + currentDirectory +"/recall.txt >> " + currentDirectory + "/recall_tmp.txt"
		subprocess.check_output(['bash','-c', cmdCat])

		cmdSed = "sed -i 's/unknown/" + soft + "/g' " + currentDirectory + "/correct_base_rate.txt"
		subprocess.check_output(['bash','-c', cmdSed])
		cmdCat = "grep " + soft + " " + currentDirectory + "/correct_base_rate.txt >> " + currentDirectory + "/correct_base_rate_tmp.txt"
		subprocess.check_output(['bash','-c', cmdCat])







def writeLatex(options, currentDirectory, errorRate, coverage, outDir, outputPDFName):
	content = r'''\documentclass{article}
	\usepackage{graphicx}

	\begin{document}
	
	\section{Recall, precision, correct bases rate}
	\begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(recall)s}
	\caption{\textbf{Recall of correctors on %(coverageToKeep)sX reads} Recall values in ordinate are computed after correction for each read experiment, using correctors in absciss. Error rate was %(errorToKeep)s)}
	\label{fig:recall}
	\end{figure}
	
	\begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(precision)s}
	\caption{\textbf{Precision of correctors on %(coverageToKeep)sX reads} Precision values in ordinate are computed after correction for each read experiment, using correctors in absciss. Error rate was %(errorToKeep)}
	\label{fig:precision}
	\end{figure}

	
	\begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(correctRate)s}
	\caption{\textbf{Correct base rate after correction on %(coverageToKeep)sX reads} Correct base rate in ordinate, computed after correction for each read experiment, using correctors in absciss. Error rate was %(errorToKeep)}
	\label{fig:correctRate}
	\end{figure}

	\section{Reads size}

	\begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(size)s}
	\caption{\textbf{Ratio of corrected over real isoform length in corrected reads} Coverage of %(coverageToKeep)sX, ratio in ordinate, corrector in absciss. Error rate was %(errorToKeep)}
	\label{fig:size}
	\end{figure}

	\section{Isoform correction}
	
	\begin{figure}[ht!]
	 \centering\includegraphics[width=0.7\textwidth]{%(isoform)s}
	 \caption{\textbf{Confusion matrix of isoforms, for coverage %(coverageToKeep)sX} Original isoforms are in absciss, corrected isoforms are in ordinate. For each pair we compute the number of corrected reads in a given isoform / original number of isoform. If no isoform was transformed during correction, the ratio is 1. Confusion matrix are presented for each correction method.}
	 \label{fig:confusionAll}
	\end{figure}'''
	
	if len(errorRate) > 1:
		content += r''' \subsection{Isoform conservation in function of error rate for exon method}
			 \begin{figure}[ht!]
		 \centering\includegraphics[width=0.7\textwidth]{%(isoformError)s}
	 \caption{\textbf{Confusion matrix of isoforms, for coverage %(coverage)sX and exon correction method.} Original isoforms are in absciss, corrected isoforms are in ordinate. For each pair we compute the number of corrected reads in a given isoform / original number of isoform. If no isoform was transformed during correction, the ratio is 1. Confusion matrix are presented for each tested error rate}
	 \label{fig:confusionError}
	\end{figure}
	 '''

	if len(coverage) > 1:
				content += r''' \subsection{Isoform conservation in function of coverage for exon method}
				\begin{figure}[ht!]
				\centering\includegraphics[width=0.7\textwidth]{%(isoformCoverage)s}
				 \caption{\textbf{Confusion matrix of isoforms, for an error rate of %(errorToKeep)s and exon correction method.} Original isoforms are in absciss, corrected isoforms are in ordinate. For each pair we compute the number of corrected reads in a given isoform / original number of isoform. If no isoform was transformed during correction, the ratio is 1. Confusion matrix are presented for each tested coverage}
				 \label{fig:confusionError}
				\end{figure}
				 '''
		

	content += r''' \section{Metrics in function of coverage for the exon-correction solution.}
	 
	 \begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(coverage_function)s}
	\caption{\textbf{Correction quality over coverage for correction by exon solution} Precision, recall and correct base ratio are colored bars. Coverages increase in the absciss.}
	\label{fig:coverage}
	\end{figure}'''

	if len(errorRate) > 1:
		content += r''' \section{Metrics in function of error rate for the exon-correction solution.}
		 \begin{figure}[ht!]
		\centering\includegraphics[width=0.8\textwidth]{%(errorrate_function)s}
		\caption{\textbf{Correction quality over error rates for correction by exon solution} Precision, recall and correct base ratio are colored bars. Error rates increase in the absciss. Computed for a coverage of %(coverageToKeep)s.}
		\label{fig:errorrate}
		\end{figure}'''
	content += r'''
	\end{document}
	'''
	print(outDir + "/" + outputPDFName +'.tex')
	with open(outDir + "/" + outputPDFName +'.tex','w') as f:
		f.write(content%options)
	proc = subprocess.Popen(['pdflatex', '-output-directory', outDir, outputPDFName + ".tex"])
	proc.communicate()



#R functions
def printConfusionMatrix(currentDirectory, corrector, confusionFile, suffix, coverageToKeep, errorToKeep):
	#~ Rcmd = "Rscript " + currentDirectory + "/plot_confusion_matrix.R " + confusionFile + " " + corrector  + " " + currentDirectory
	Rcmd = "Rscript " + currentDirectory + "/plot_all_confusion_matrix.R " + currentDirectory + "/all_confusion_matrix.txt " + currentDirectory + " " + str(coverageToKeep) + " " + str(errorToKeep)
	subprocessLauncher(Rcmd)


def printConfusionMatrixFunctionOf(currentDirectory, coverageToKeep, errorToKeep):
	if errorToKeep is not None:
		Rcmd = "Rscript " + currentDirectory + "/plot_confusion_matrix_function_error.R " + currentDirectory + "/all_confusion_matrix.txt " + currentDirectory + " "  + str(coverageToKeep)
		subprocessLauncher(Rcmd)
	if coverageToKeep is not None:
		Rcmd = "Rscript " + currentDirectory + "/plot_confusion_matrix_function_coverage.R " + currentDirectory + "/all_confusion_matrix.txt " + currentDirectory + " "  + str(errorToKeep) 
		subprocessLauncher(Rcmd)
	

def printMetrics(currentDirectory, cov, error):
	cmdR = "Rscript " + currentDirectory + "/plot_recall.R " + currentDirectory + "/recall_cov_" + str(cov) + "_err_" + str(error) + ".txt " + currentDirectory
	subprocessLauncher(cmdR)
	cmdR = "Rscript " + currentDirectory + "/plot_precision.R " + currentDirectory + "/precision_cov_" + str(cov) + "_err_" + str(error) +".txt " + currentDirectory
	subprocessLauncher(cmdR)
	cmdR = "Rscript " + currentDirectory + "/plot_correct_base_rate.R " + currentDirectory + "/correct_base_rate_cov_" + str(cov)  + "_err_" + str(error) +".txt " + currentDirectory
	subprocessLauncher(cmdR)
	cmdR = "Rscript " + currentDirectory + "/plot_size.R " + currentDirectory + "/sizes_reads_cov_" + str(cov) + "_err_" + str(error) +".txt " + currentDirectory
	subprocessLauncher(cmdR)
	cmdR = "Rscript " + currentDirectory + "/plot_all_metrics_coverage.R " + currentDirectory + "/all_recall_precision.txt " + currentDirectory
	subprocessLauncher(cmdR)


def printMetricErrorRates(currentDirectory, cov):
	cmdR = "Rscript " + currentDirectory + "/plot_all_metrics_errorrate.R " + currentDirectory + "/all_recall_precision.txt " + str(cov) + " " + currentDirectory
	subprocessLauncher(cmdR)

def computeResultsIsoforms(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts, cov, covToPrint, allCoverages, errorRate, errorToPrint, allErrorRates, outDir="/home/marchet/detection-consensus-isoform/results"):
	msa(suffix, correc)
	correcIsoformTypesToCounts, correcIsoformTypesToSeq = makeCorrectedHeadersList(outDir, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts)
	confusionFile, isCorrect = computeRatioIsoforms(refIsoformTypesToCounts, correcIsoformTypesToCounts, currentDirectory, suffix, correc, cov, errorRate)
	if cov == covToPrint and errorRate == errorToPrint:
		printConfusionMatrix(currentDirectory, correc, confusionFile, suffix, cov, errorRate)
	if len(allCoverages) > 1:
		printConfusionMatrixFunctionOf(currentDirectory, covToPrint, errorToPrint)
	if len(allErrorRates) > 1:
		printConfusionMatrixFunctionOf(currentDirectory, covToPrint, errorToPrint)
	return isCorrect





def checkReadFiles(readfiles):
	if readfiles is None:
		return True
	allFilesAreOK = True
	#~ for file in readfiles:
	if not os.path.isfile(readfiles):
		print("[ERROR] File \""+readfiles+"\" does not exist.")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more read files do not exist.")


def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))

	print(currentDirectory)
	
	covSR = 1
	#~ skipped = [40,50,100]
	skipped = [100]
	#~ abund = [50,75,90,10]
	#~ abund = [50,75,90]
	#~ errorRate = [13,5]
	errorRate = [5,9,13,15]
	errorRateToKeep = 13

	abund = [50,75,90]
	#~ skipped = [100]
	skippedS = [str(r) for r in skipped]
	abundS = [str(r) for r in abund]
	EStype = "ES"
	#~ EStype = "MES"
	#~ EStype = "alt"

	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	parser = argparse.ArgumentParser()
	parser.add_argument('-output', nargs='?', type=str, action="store", dest="outputDirPath", help="Name for output directory", default=None)
	#~ parser.add_argument('-corrector', type=str, action="store", dest="correctors", help="A particular corrector to be used", default=None)
	parser.add_argument('-coverage', nargs='?', type=int, action="store", dest="coverage", help="Coverage for LR simulation (default 20)", default=None)
	parser.add_argument('-pdf', nargs='?', type=str, action="store", dest="pdf", help="Name for pdf outfile", default="expe")
	# get options for this run
	args = parser.parse_args()
	outputDirPath = args.outputDirPath
	#~ covLR = args.covLR

	#order to keep for plots
	#LoRMA",
     #LoRDEC",
     # "Proovread",
     #"MECAT",
     #"PBDagCon",
     #"daccord",
     #"msa_exon",
     #"msa_isoform",
     #msa_both
     # "colorMap"
	#~ correctors = ["msa_exon","msa_isoform"]
	#~ correctors = ["msa_isoform"]
	correctors = ["msa_exon", "msa_isoform", "msa_both"]
	#~ correctors = ["msa_both"]
	#~ correctors = ["msa_sparc"]
	if args.coverage is not None:
		coverage = [args.coverage]
	else:
		#~ coverage = [15]
		#~ coverage = [5, 10, 15, 20, 30, 50, 100]
		#~ coverage = [10, 15, 20, 30, 50]
		coverage = [10, 15, 20]
	covForFigs = 15
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
	else:
		outputDirPath = currentDirectory
	outputPDFName = args.pdf 

	cmdFile = '''echo "value metric coverage error" > ''' + currentDirectory + '''/all_recall_precision.txt'''
	subprocess.check_output(['bash','-c', cmdFile])
	cmdFile = '''echo "reference correction ratio soft coverage error" > ''' + currentDirectory + '''/all_confusion_matrix.txt'''
	subprocess.check_output(['bash','-c', cmdFile])
	


	simulateReads(skipped, abund, coverage, EStype, currentDirectory, errorRate)
	for error in errorRate:
		for covLR in coverage:
			#~ input(covLR)
			print(covLR)
			outSize = open(currentDirectory + "/sizes_reads_cov_"+ str(covLR)+ "_err_" +str(error)+".txt", 'w')
			outSize.write("soft size skipped abund\n")
			for correc in correctors:
				for skippedExon in skippedS:
					for abundanceMajor in abundS:
						listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq = makeReferenceHeadersList(currentDirectory, str(skippedExon), str(abundanceMajor), str(covLR), str(error))
						suffix = "_size_" + str(skippedExon) + "_abund_" + str(abundanceMajor) + "_cov_" + str(covLR) + "_err_" + str(error)
						isCorrect = computeResultsIsoforms(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts, covLR, covForFigs, coverage, error,errorRateToKeep, errorRate)
						if isCorrect: # all reads were corrected to the right isoform
							alignOnRefMsa(correc, skippedExon, abundanceMajor, currentDirectory, "/home/marchet/detection-consensus-isoform/results", covLR, error)
							computeResultsRecallPrecision(correc, skippedExon, abundanceMajor, currentDirectory, correc, refIsoformTypesToSeq, outSize, covLR, error)
			cmdMv = "mv " + currentDirectory + "/recall_tmp.txt " + currentDirectory + "/recall_cov_"+ str(covLR) + "_err_" + str(error) + ".txt"
			subprocess.check_output(['bash','-c', cmdMv])
			cmdMv = "mv " + currentDirectory + "/precision_tmp.txt " + currentDirectory + "/precision_cov_"+ str(covLR)+ "_err_" + str(error)  +".txt"
			subprocess.check_output(['bash','-c', cmdMv])
			outSize.close()
			cmdMv = "mv " + currentDirectory + "/correct_base_rate_tmp.txt " + currentDirectory + "/correct_base_rate_cov_"+ str(covLR) + "_err_" + str(error) +".txt"
			subprocess.check_output(['bash','-c', cmdMv])
			cmdAwk = '''awk '{if($1=="msa_exon") {print$2, "recall",''' + str(covLR) + ''',''' +str(error)  + '''}}' ''' + currentDirectory + '''/recall_cov_'''+ str(covLR) + '''_err_''' + str(error) +'''.txt >>''' + currentDirectory +  '''/all_recall_precision.txt'''
			print(cmdAwk)
			subprocess.check_output(['bash','-c', cmdAwk])
			cmdAwk = '''awk '{if($1=="msa_exon") {print$2, "precision",''' + str(covLR) + ''',''' +str(error) + '''}}' ''' + currentDirectory + '''/precision_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_recall_precision.txt'''
			print(cmdAwk)
			subprocess.check_output(['bash','-c', cmdAwk])
			cmdAwk = '''awk '{if($1=="msa_exon") {print$2, "correct",''' + str(covLR) + ''',''' +str(error) +'''}}' ''' + currentDirectory + '''/correct_base_rate_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_recall_precision.txt'''
			print(cmdAwk)
			subprocess.check_output(['bash','-c', cmdAwk])
			if error == errorRateToKeep and covLR == covForFigs:
				printMetrics(currentDirectory, covForFigs, errorRateToKeep)
	printMetricErrorRates(currentDirectory, covForFigs)
	dictLatex = {"coverage":str(covLR), "recall": currentDirectory + "/recall.png", "precision": currentDirectory + "/precision.png", "correctRate": currentDirectory + "/correct_base_rate.png", "size":  currentDirectory + "/size.png", "coverage_function": currentDirectory + "/metrics_function_coverage.png", "errorrate_function" : currentDirectory + "/metrics_function_errorrate.png", "coverageToKeep": covForFigs, "errorToKeep": errorRateToKeep, "isoform": currentDirectory + "/all_confusion_matrix.png", "isoformError": currentDirectory+"/confusion_matrix_function_error.png", "isoformCoverage": currentDirectory+"/confusion_matrix_function_coverage.png"}
	writeLatex(dictLatex, currentDirectory, errorRate, coverage, outputDirPath, outputPDFName)

if __name__ == '__main__':
	main()
