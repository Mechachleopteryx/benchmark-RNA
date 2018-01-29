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
def simulateReads(covSR, covLR, skipped, abund, EStype, currentDirectory):
	for sizeSkipped in skipped:
		for relAbund in abund:
			suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
			# simulation
			if EStype == "ES":
				cmdSimul = currentDirectory + "/ES_simulation " + str(sizeSkipped) + " " + str(relAbund) + " " + str(suffix) + " " +  str(covSR) + " " + str(covLR)
			elif EStype == "MES":
				cmdSimul =  currentDirectory + "/MES_simulation " + str(sizeSkipped) + " " + str(relAbund) + " " + str(suffix) +  " " + str(covLR)
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
	isCorrect = True
	for ref in counts:
		for ref2 in counts:
			if ref2 in counts[ref].keys():
				ratio = counts[ref][ref2] * 1.0 / len(refIsoformTypesToCounts[ref]) if len(refIsoformTypesToCounts[ref]) != 0 else 0
				if ratio != 1:
					isCorrect = False
				outConf.write(ref + " " + ref2 + " " + str(ratio) + "\n")
			#~ else:
				#~ print("this", ref2)
				#~ outConf.write(ref + " " + ref2 + " 0\n")
				#~ isCorrect = False
	outConf.close()
	return confusionName, isCorrect



	


# check results by aligning on the reference sequences
#~ def alignOnRef(listCorrectors, skipped, abund, outDir="/home/marchet/detection-consensus-isoform/results"):
	#~ for soft in listCorrectors:
		#~ for sizeSkipped in skipped:
			#~ for relAbund in abund:
				#~ suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
				#~ samFile = open("results" + soft + suffix + ".sam", 'w')
				#~ cmdAlign = "/home/marchet/bin/Complete-Striped-Smith-Waterman-Library/src/ssw_test refSequences" + suffix + ".fa " + outDir + "/corrected_by_" + soft + suffix + ".fa  -c -s"
				#~ p = subprocessLauncher(cmdAlign, samFile)




def alignOnRefMsa(soft, skipped, abund, currentDirectory, resultDirectory):
	suffix = "_size_" + str(skipped) + "_abund_" + str(abund)
	listFileNames = getFiles(resultDirectory, "corrected_by_MSA*.fa")
	for fileC in listFileNames:
		isoform = getCorrectedHeaders(resultDirectory + "/" + fileC)[0].split("_")[0][1:]
		print(isoform, "**********************isoform")
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



def getExpectedLength(currentDirectory, suffix):
	expectedLengths = {"exclusion":0, "inclusion":0}
	refFile = open(currentDirectory + "/refSequences" + suffix + ".fa", 'r')
	lines = refFile.readlines()
	for l in lines:
		if ">" in l:
			targetType = l[1:-1]
		else:
			expectedLengths[targetType] = len(l) - 1
	return expectedLengths


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



def computeResultsRecallPrecision(corrector, skipped, abund, currentDirectory, soft, refIsoformTypesToSeq):
	print("********************************************************")
	suffix = "_size_" + str(skipped) + "_abund_" + str(abund)
	expectedLengths = getExpectedLength(currentDirectory, suffix)
	meanSizes = {}
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", soft)
	for isofType in refIsoformTypesToSeq:
		print(isofType)
		start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries = readSam(soft, suffix, isofType, currentDirectory)
		meanSizes[isofType] = {"realSize" : [], "alignedSize" : []}
		for querySeq, aln in blockResults.items():
		# meanSizes = {"exclusion":{"realSize":[], "alignedSize":[]}, "inclusion":{"realSize":[], "alignedSize":[]}}
			#~ meanSizes[isofType] = {}
			meanSizes[isofType]["realSize"].append(lenResults[querySeq][isofType][0])
			meanSizes[isofType ]["alignedSize"].append(lenResults[querySeq][isofType][1])
		print(currentDirectory + "/corrected_reads_by_" + soft + "_" + isofType + suffix + ".fa #####################################################################################")
		cmdHead = "head -2 " + currentDirectory + "/corrected_reads_by_" + soft + "_" + isofType + suffix + ".fa > " + currentDirectory + "/corrected.fa"
		subprocess.check_output(['bash','-c', cmdHead])
		cmdHead = "head -2 " + currentDirectory + "/uncorrected_reads_"  + isofType + suffix + ".fa > " + currentDirectory + "/uncorrected.fa"
		subprocess.check_output(['bash','-c', cmdHead])
		cmdHead = "head -2 " + currentDirectory + "/perfect_reads_" + isofType + suffix + ".fa > " + currentDirectory + "/perfect.fa"
		subprocess.check_output(['bash','-c', cmdHead])
		cmdBench = "python3 " + currentDirectory + "/benchmark-long-read-correction/benchmark.py -c " + currentDirectory + "/corrected.fa -u " + currentDirectory + "/uncorrected.fa -r " + currentDirectory + "/perfect.fa -o " + currentDirectory
		subprocessLauncher(cmdBench)
		#~ cmdMv = "mv " + currentDirectory + "/outPositions.txt " + currentDirectory + "/outPositionsExclusion_" + soft +  suffix + ".fa"
		#~ subprocess.check_output(['bash','-c', cmdMv])

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

	

	#~ if start is not None:
	#~ nbIncl = getFileReadNumber("perfect_reads_inclusion" + suffix + ".fa")
	#~ nbExcl = getFileReadNumber("perfect_reads_exclusion" + suffix + ".fa")
	#~ meanSizes, countIncl, countExcl, countAssigned, countInverted = getIsoform(blockResults, lenResults, suffix, queries, nbIncl, nbExcl, soft)
	#~ readNumber = len(blockResults.keys())
	meanReadsSize = round(sum(readsSize)*1.0/len(readsSize),2) if len(readsSize) > 0 else 0
	meanInclusionCorrectedSize = round(sum(meanSizes["inclusion"]["alignedSize"])*1.0/len(meanSizes["inclusion"]["alignedSize"]),2) if len(meanSizes["inclusion"]["alignedSize"]) > 0 else 0
	meanExclusionCorrectedSize = round(sum(meanSizes["exclusion"]["alignedSize"])*1.0/len(meanSizes["exclusion"]["alignedSize"]),2) if len(meanSizes["exclusion"]["alignedSize"]) > 0 else 0

	ratioLenI = round(meanInclusionCorrectedSize*100/expectedLengths["inclusion"],2)
	ratioLenE = round(meanExclusionCorrectedSize*100/expectedLengths["exclusion"],2)

	print("#############################")
	#~ print(soft, "correction for inclusion size", sizeSkipped, ", ratio:", relAbund, "/", 100-relAbund)
	#~ print("Corrected inclusion reads in output:", round(countIncl*100.0/(countIncl+countExcl),2) if countIncl+countExcl != 0 else 0, "%")
	#~ print("Corrected exclusion reads in output:", round(countExcl*100.0/(countExcl+countIncl),2) if countIncl+countExcl != 0 else 0, "%")
	#~ print("Mean length corrected reads:", meanReadsSize, ", mean length of inclusion reads:", meanInclusionCorrectedSize, "(", ratioLenI, "% of real size ), mean length of exclusion reads:", meanExclusionCorrectedSize, "(", ratioLenE, "% of real size )")
	#~ outSize.write(soft + " " + str(ratioLenI) + " " +  str(sizeSkipped) + " "+ str(relAbund) + "\n")
	#~ outSize.write(soft + " " + str(ratioLenE) + " " +  str(sizeSkipped) + " "+ str(relAbund) +"\n")
	#~ correctedToIncl = round(countIncl*100.0/(countIncl+countExcl),2) if countIncl+countExcl != 0 else 0
	#~ if not correctedToIncl == 0 :
		#~ out.write(soft + " " + str(sizeSkipped) + " " + str(relAbund) + " " + str(correctedToIncl) +"\n")
		#~ out2.write(soft + " inclusion " +  str(relAbund)  + " " + str(ratioLenI) + " " + str(sizeSkipped) +"\n")
		#~ out2.write(soft + " exclusion " +  str(relAbund)  + " " + str(ratioLenE) + " " + str(sizeSkipped) +"\n")

	### launch corrector benchmark"
	#~ if "msa" not in soft :
		#~ if countExcl > 0:
			#~ cmdBench = "python3 ./benchmark-long-read-correction/benchmark.py -c corrected_reads_exclusion" + suffix + ".fa -u uncorrected_reads_exclusion" + suffix + ".fa -r perfect_reads_exclusion" + suffix + ".fa"
			#~ print("Exclusion")
			#~ subprocessLauncher(cmdBench)
			#~ cmdMv = "mv outPositions.txt outPositionsExclusion_" + soft + suffix + ".fa"
			#~ subprocess.check_output(['bash','-c', cmdMv])

			#~ cmdSed = "sed -i 's/unknown/" + soft + "/g' precision.txt"
			#~ subprocess.check_output(['bash','-c', cmdSed])
			#~ cmdCat = "grep " + soft + " precision.txt >> precision_tmp.txt"
			#~ subprocess.check_output(['bash','-c', cmdCat])
			
			#~ cmdSed = "sed -i 's/unknown/" + soft + "/g' recall.txt"
			#~ subprocess.check_output(['bash','-c', cmdSed])
			#~ cmdCat = "grep " + soft + " recall.txt >> recall_tmp.txt"
			#~ subprocess.check_output(['bash','-c', cmdCat])
			
			#~ cmdSed = "sed -i 's/unknown/" + soft + "/g' correct_base_rate.txt"
			#~ subprocess.check_output(['bash','-c', cmdSed])
			#~ cmdCat = "grep " + soft + " correct_base_rate.txt >> correct_base_rate_tmp.txt"
			#~ subprocess.check_output(['bash','-c', cmdCat])
			
		#~ if countIncl > 0:
			#~ cmdBench = "python3 ./benchmark-long-read-correction/benchmark.py -c corrected_reads_inclusion" + suffix + ".fa -u uncorrected_reads_inclusion" + suffix + ".fa -r perfect_reads_inclusion" + suffix + ".fa"
			#~ print("Inclusion")
			#~ subprocessLauncher(cmdBench)
			#~ cmdMv = "mv outPositions.txt outPositionsInclusion_" + soft + suffix + ".fa"
			#~ subprocess.check_output(['bash','-c', cmdMv])

			#~ cmdSed = "sed -i 's/unknown/" + soft + "/g' precision.txt"
			#~ subprocess.check_output(['bash','-c', cmdSed])
			#~ cmdCat = "grep " + soft + " precision.txt >> precision_tmp.txt"
			#~ subprocess.check_output(['bash','-c', cmdCat])
			
			#~ cmdSed = "sed -i 's/unknown/" + soft + "/g' recall.txt"
			#~ subprocess.check_output(['bash','-c', cmdSed])
			#~ cmdCat = "grep " + soft + " recall.txt >> recall_tmp.txt"
			#~ subprocess.check_output(['bash','-c', cmdCat])

			#~ cmdSed = "sed -i 's/unknown/" + soft + "/g' correct_base_rate.txt"
			#~ subprocess.check_output(['bash','-c', cmdSed])
			#~ cmdCat = "grep " + soft + " correct_base_rate.txt >> correct_base_rate_tmp.txt"
			#~ subprocess.check_output(['bash','-c', cmdCat])
	#~ else:
		#~ print("Exclusion")
		#~ if countExcl > 0:
	




#R functions
def printConfusionMatrix(currentDirectory, corrector, confusionFile, suffix):
	Rcmd = "Rscript " + currentDirectory + "/matrice_confusion.R " + confusionFile + " " + corrector  + suffix + " " + currentDirectory
	subprocessLauncher(Rcmd)


def computeResultsIsoforms(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts, outDir="/home/marchet/detection-consensus-isoform/results"):
	msa(suffix, correc)
	correcIsoformTypesToCounts, correcIsoformTypesToSeq = makeCorrectedHeadersList(outDir, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts)
	confusionFile, isCorrect = computeRatioIsoforms(refIsoformTypesToCounts, correcIsoformTypesToCounts, currentDirectory, suffix)
	printConfusionMatrix(currentDirectory, correc, confusionFile, suffix)
	return isCorrect


def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))
	covSR = 1
	#~ skipped = [50,100]
	#~ abund = [50,75,90,10]
	abund = [50]
	skipped = [100]
	skippedS = [str(r) for r in skipped]
	abundS = [str(r) for r in abund]
	EStype = "ES"
	#~ EStype = "MES"

	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	parser = argparse.ArgumentParser()
	parser.add_argument('-output', nargs='?', type=str, action="store", dest="outputDirPath", help="Name for output directory", default=None)
	#~ parser.add_argument('-corrector', type=str, action="store", dest="correctors", help="A particular corrector to be used", default=None)
	parser.add_argument('-coverage', nargs='?', type=int, action="store", dest="covLR", help="Coverage for LR simulation (default 10)", default=20)
	# get options for this run
	args = parser.parse_args()
	outputDirPath = args.outputDirPath
	covLR = args.covLR

	correctors = ["msa_isoform", "msa_exon"]
	#~ correctors = ["msa_isoform"]
	#~ correctors = ["msa_exon"]
	#~ correctors = ["msa_sparc"]
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
				
	simulateReads(covSR, covLR, skipped, abund, EStype, currentDirectory)
	for correc in correctors:
		for skippedExon in skippedS:
			for abundanceMajor in abundS:
				listFilesPerfect, refIsoformTypesToCounts, refIsoformTypesToSeq = makeReferenceHeadersList(currentDirectory, str(skippedExon), str(abundanceMajor))
				suffix = "_size_" + str(skippedExon) + "_abund_" + str(abundanceMajor)
				isCorrect = computeResultsIsoforms(correc, currentDirectory, skippedExon, abundanceMajor, suffix, refIsoformTypesToCounts)
				if isCorrect: # all reads were corrected to the right isoform
					alignOnRefMsa(correc, skippedExon, abundanceMajor, currentDirectory, "/home/marchet/detection-consensus-isoform/results")
					computeResultsRecallPrecision(correc, skippedExon, abundanceMajor, currentDirectory, correc, refIsoformTypesToSeq)
	cmdMv = "mv " + currentDirectory + "/recall_tmp.txt " + currentDirectory + "/recall.txt"
	subprocess.check_output(['bash','-c', cmdMv])
	cmdMv = "mv " + currentDirectory + "/precision_tmp.txt " + currentDirectory + "/precision.txt"
	subprocess.check_output(['bash','-c', cmdMv])

	cmdMv = "mv " + currentDirectory + "/correct_base_rate_tmp.txt " + currentDirectory + "/correct_base_rate.txt"
	subprocess.check_output(['bash','-c', cmdMv])

if __name__ == '__main__':
	main()
