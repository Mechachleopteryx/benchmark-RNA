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
from utils import *




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
	listFileNames = getFiles(resultDirectory, "corrected_by_MSA*.fa")#todo change here MSA => specific
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







def readSam(soft, suffix, pathSam, currentDirectory):
	blockResults = dict()
	lenResults = dict()
	#~ pathSam = currentDirectory + "/results" +isoformType + soft + suffix + ".sam"
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
							
		pathSam = currentDirectory + "/results" + isofType + soft + suffix + ".sam"

		start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries = readSam(soft, suffix, pathSam, currentDirectory)
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


def benchMsa(errorRate, coverage, correctors, skippedS, abundS, currentDirectory, covForFigs, errorRateToKeep, outputDirPath, outputPDFName):
	for error in errorRate:
		for covLR in coverage:
			outSize = open(currentDirectory + "/sizes_reads_cov_"+ str(covLR)+ "_err_" +str(error)+".txt", 'w')
			outSize.write("soft size skipped abund\n")
			for correc in correctors:
				for skippedExon in skippedS:
					for abundanceMajor in abundS:
						listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq = makeReferenceHeadersList(currentDirectory, str(skippedExon), str(abundanceMajor), str(covLR), str(error))

						if "msa" in correc:
							#~ listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq = makeReferenceHeadersList(currentDirectory, str(skippedExon), str(abundanceMajor), str(covLR), str(error))
							suffix = "_size_" + str(skippedExon) + "_abund_" + str(abundanceMajor) + "_cov_" + str(covLR) + "_err_" + str(error)
							isCorrect = computeResultsIsoforms(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts, covLR, covForFigs, coverage, error,errorRateToKeep, errorRate)
							if isCorrect: # all reads were corrected to the right isoform
								alignOnRefMsa(correc, skippedExon, abundanceMajor, currentDirectory, "/home/marchet/detection-consensus-isoform/results", covLR, error)
								computeResultsRecallPrecision(correc, skippedExon, abundanceMajor, currentDirectory, correc, refIsoformTypesToSeq, outSize, covLR, error)
						else:
							launchCorrector(currentDirectory, correc , skippedExon, abundanceMajor, covLR, error)
							alignOnRef(correc,skippedExon, abundanceMajor, currentDirectory, covLR, error)
							computeResults(correc,skippedExon, abundanceMajor, currentDirectory, covLR, error, refIsoformTypesToSeq)

			cmdMv = "mv " + currentDirectory + "/recall_tmp.txt " + currentDirectory + "/recall_cov_"+ str(covLR) + "_err_" + str(error) + ".txt"
			subprocess.check_output(['bash','-c', cmdMv])
			cmdMv = "mv " + currentDirectory + "/precision_tmp.txt " + currentDirectory + "/precision_cov_"+ str(covLR)+ "_err_" + str(error)  +".txt"
			subprocess.check_output(['bash','-c', cmdMv])
			outSize.close()
			cmdMv = "mv " + currentDirectory + "/correct_base_rate_tmp.txt " + currentDirectory + "/correct_base_rate_cov_"+ str(covLR) + "_err_" + str(error) +".txt"
			subprocess.check_output(['bash','-c', cmdMv])

			
			cmdAwk = '''awk '{print $0,''' + str(covLR) + ''',''' +str(error) + '''}' ''' + currentDirectory + '''/precision_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_precisions_softs.txt'''
			subprocess.check_output(['bash','-c', cmdAwk])
			cmdAwk = '''awk '{print $0,''' + str(covLR) + ''',''' +str(error) + '''}' ''' + currentDirectory + '''/recall_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_recalls_softs.txt'''
			subprocess.check_output(['bash','-c', cmdAwk])
			cmdAwk = '''awk '{print $0,''' + str(covLR) + ''',''' +str(error) + '''}' ''' + currentDirectory + '''/correct_base_rate_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_correctRates_softs.txt'''
			subprocess.check_output(['bash','-c', cmdAwk])

			
			cmdAwk = '''awk '{if($1=="msa_exon") {print$2, "recall",''' + str(covLR) + ''',''' +str(error)  + '''}}' ''' + currentDirectory + '''/recall_cov_'''+ str(covLR) + '''_err_''' + str(error) +'''.txt >>''' + currentDirectory +  '''/all_recall_precision.txt'''
			subprocess.check_output(['bash','-c', cmdAwk])
			cmdAwk = '''awk '{if($1=="msa_exon") {print$2, "precision",''' + str(covLR) + ''',''' +str(error) + '''}}' ''' + currentDirectory + '''/precision_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_recall_precision.txt'''
			subprocess.check_output(['bash','-c', cmdAwk])
			cmdAwk = '''awk '{if($1=="msa_exon") {print$2, "correct",''' + str(covLR) + ''',''' +str(error) +'''}}' ''' + currentDirectory + '''/correct_base_rate_cov_'''+ str(covLR) + '''_err_''' + str(error)+ '''.txt >>''' + currentDirectory +  '''/all_recall_precision.txt'''
			subprocess.check_output(['bash','-c', cmdAwk])
			if error == errorRateToKeep and covLR == covForFigs:
				printMetrics(currentDirectory, covForFigs, errorRateToKeep)
	printGlobalMetrics(currentDirectory)
	printMetricErrorRates(currentDirectory, covForFigs)
	dictLatex = {"coverage":str(covLR), "recall": currentDirectory + "/all_recall_softs.png", "precision": currentDirectory + "/all_precision_softs.png", "correctRate": currentDirectory + "/all_correctRate_softs.png", "size":  currentDirectory + "/size.png", "coverage_function": currentDirectory + "/metrics_function_coverage.png", "errorrate_function" : currentDirectory + "/metrics_function_errorrate.png", "coverageToKeep": covForFigs, "errorToKeep": errorRateToKeep, "isoform": currentDirectory + "/all_confusion_matrix.png", "isoformError": currentDirectory+"/confusion_matrix_function_error.png", "isoformCoverage": currentDirectory+"/confusion_matrix_function_coverage.png", "precisionSoft": currentDirectory+"/plot_all_precision_softs.png", "recallSoft": currentDirectory+"/plot_all_recall_softs.png"}
	writeLatex(dictLatex, currentDirectory, errorRate, coverage, outputDirPath, outputPDFName)
