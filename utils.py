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



######### utils for warnings, user messages, errors ######### 

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

######### utils for subprocess ######### 



# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
	return p


######### utils for sequence files ######### 


# find files with a regex
def getFiles(pathToFiles, name): #for instance name can be "*.txt"
	os.chdir(pathToFiles)
	listFiles = []
	for files in glob.glob(name):
		listFiles.append(files)
	return listFiles

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



######### utils for read simulation ######### 
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




######### utils for latex ######### 


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
	\caption{\textbf{Precision of correctors on %(coverageToKeep)sX reads} Precision values in ordinate are computed after correction for each read experiment, using correctors in absciss. Error rate was %(errorToKeep)s)}
	\label{fig:precision}
	\end{figure}
	
	\begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(correctRate)s}
	\caption{\textbf{Correct base rate after correction %(coverageToKeep)sX reads} Precision values in ordinate are computed after correction for each read experiment, using correctors in absciss. Error rate was %(errorToKeep)s)}
	\label{fig:correctRate}
	\end{figure}


	\section{Reads size}
	\begin{figure}[ht!]
	\centering\includegraphics[width=0.8\textwidth]{%(size)s}
	\caption{\textbf{Ratio of corrected over real read size after correction} Coverage of %(coverageToKeep)sX , ratio in ordinate, correctors in absciss. Error rate was %(errorToKeep)s}
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
	\end{figure}'''

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
	\caption{\textbf{Correction quality over coverage for correction by exon solution} Mean precision, recall and correct base ratio are colored bars. Coverages increase in the absciss.}
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
		\end{document} '''

	print(outDir + "/" + outputPDFName +'.tex')
	with open(outDir + "/" + outputPDFName +'.tex','w') as f:
		f.write(content%options)
	proc = subprocess.Popen(['pdflatex', '-output-directory', outDir, outputPDFName + ".tex"])
	proc.communicate()





######### utils for R functions ######### 

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


def printGlobalMetrics(currentDirectory):
	cmdR = "Rscript " + currentDirectory + "/plot_all_precisions_softs.R " + currentDirectory + "/all_precisions_softs.txt " + currentDirectory
	subprocessLauncher(cmdR)
	cmdR = "Rscript " + currentDirectory + "/plot_all_recalls_softs.R " + currentDirectory + "/all_recalls_softs.txt " + currentDirectory
	subprocessLauncher(cmdR)
	cmdR = "Rscript " + currentDirectory + "/plot_all_correctRate_softs.R " + currentDirectory + "/all_correctRates_softs.txt " + currentDirectory
	subprocessLauncher(cmdR)


def printMetricErrorRates(currentDirectory, cov):
	cmdR = "Rscript " + currentDirectory + "/plot_all_metrics_errorrate.R " + currentDirectory + "/all_recall_precision.txt " + str(cov) + " " + currentDirectory
	subprocessLauncher(cmdR)
