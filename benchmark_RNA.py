#!/usr/bin/env python3
import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import re
import copy

import seaborn as sn
import pandas as pd
import numpy as np
#~ import matplotlib.pyplot as plt
#~ from sklearn.metrics import confusion_matrix

def readfasta(infile):
    labels = []
    sequences = []

    curlabel = None
    cursequence = ""

    def updatelists():
        if len(cursequence) is not 0:
            sequences.append(cursequence)
            if curlabel is not None:
                labels.append(curlabel)
            else:
                labels.append('seq'+str(len(sequences)))

    for line in infile:
        if line[0] == ">":
            updatelists()
            cursequence = ""
            curlabel = line[1:].strip()
        else:
            cursequence += line.strip()

    updatelists()
    return list(zip(labels, sequences))



# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	#~ p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
	return p

def checkWrittenFiles(files):
	allFilesAreOK = True
	if not os.path.isfile(files):
		print("[ERROR] There was a problem writing \"" + files + "\".")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more files could not be written.")



def simulation(sizeSkipped, relAbund, suffix, covSR = 100, covLR = 10):
	cmdSimul = "./ES_simulation " + sizeSkipped + " " + relAbund + " " + suffix + " " +  str(covSR) + " " + str(covLR)
	cmdSimul = subprocessLauncher(cmdSimul)


def msa(suffix, msaType):
	if msaType == "msa_isoform":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py simulatedLR" + suffix + ".fa"
	elif msaType == "msa_exon":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py simulatedLR" + suffix + ".fa -c exon"
	elif msaType == "msa_sparc":
		cmdMSA = "/home/marchet/detection-consensus-isoform/analyze_MSAv2.py simulatedLR" + suffix + ".fa -s True"
	p = subprocessLauncher(cmdMSA)
	cmdmv = "cp /home/marchet/detection-consensus-isoform/results/corrected_by_MSA.fa corrected_by_" + msaType + suffix + ".fa"
	try:
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	cmdmv = "cp  /home/marchet/detection-consensus-isoform/results/corrected_by_MSA0.fa corrected_by_" + msaType + "_c0" + suffix + ".fa"
	try:
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	cmdmv = "cp  /home/marchet/detection-consensus-isoform/results/corrected_by_MSA1.fa corrected_by_" + msaType + "_c1" + suffix + ".fa"
	try:
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass


		
def lordec(suffix):
	cmdLordec = "lordec-correct -i simulatedLR" + suffix + ".fa -2 simulatedSR" + suffix + ".fa -o corrected_by_LoRDEC" + suffix + ".fa -k 21 -s 2 -T 4"
	p = subprocessLauncher(cmdLordec)


def colormap(suffix):
	cmdColorMap = "runCorr.sh simulatedLR"  + suffix + ".fa simulatedSR"  + suffix + ".fa colorMapDir pre 4"
	p = subprocessLauncher(cmdColorMap)
	cmdmv = "mv colorMapDir/pre_sp.fasta corrected_by_colorMap" + suffix + ".fa"
	subprocess.check_output(['bash','-c', cmdmv])


def lorma(suffix):
	cmdcp = "cp simulatedLR" + suffix + ".fa copy.fa"
	subprocess.check_output(['bash','-c', cmdcp])
	cmdLorma = "lorma.sh -threads 4 copy.fa -s"
	p = subprocessLauncher(cmdLorma)
	try:
		cmdmv = "mv trim.fasta corrected_by_LoRMA" + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	try:
		cmdmv = "mv discarded.fasta discarded_by_LoRMA" + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	try:
		cmdmv = "mv final.fasta corrected_by_LoRMA" + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		cmdnull = "> corrected_by_LoRMA" + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdnull])

def mecat(suffix):
	cmdMecat = "mecat2pw -j 0 -d simulatedLR" + suffix + ".fa  -o candidates" + suffix + ".txt  -w . -t 4 -x 1"
	p = subprocessLauncher(cmdMecat)
	cmdMecat = "mecat2cns -i 0 -t 4 -x 1 candidates" + suffix + ".txt simulatedLR" + suffix + ".fa corrected_by_MECAT" + suffix + ".fa" #-l 500
	p = subprocessLauncher(cmdMecat)


def dazzdb_daligner(suffix):
	cmdHeaders = "headersPB.pl simulatedLR" + suffix + ".fa"
	fastaPacbioHeader = open("simulatedLR" + suffix + ".fasta", 'w') #dazzdb looks for files ending exactly by .fasta (no .fa)
	p = subprocessLauncher(cmdHeaders, fastaPacbioHeader)
	fastaPacbioHeader.close()
	#launch dazzdb
	cmdDazzDB = "fasta2DB -v db_simulatedLR" + suffix  + " simulatedLR" + suffix + ".fasta"
	p = subprocessLauncher(cmdDazzDB)
	#daligner has to be launched from its directory
	#daligner binary copied from pddagcon submodules
	cmdDaligner = "daligner db_simulatedLR" + suffix + " db_simulatedLR" + suffix
	p = subprocessLauncher(cmdDaligner)
	cmdsort = "LAsort db_simulatedLR" + suffix + ".db_simulatedLR" + suffix + "*"
	subprocess.check_output(['bash','-c', cmdsort])
	cmdmerge = "LAmerge db_simulatedLR" + suffix + ".las db_simulatedLR" + suffix + ".db_simulatedLR" + suffix + ".*S.las"
	subprocess.check_output(['bash','-c', cmdmerge])



def pbdagcon(suffix):
	outDagcon = open("corrected_by_PBDagCon" + suffix + ".fa", 'w')
	cmdDagcon = "dazcon -ox -j 4 -s db_simulatedLR" + suffix + ".db -a db_simulatedLR" + suffix + ".las"
	p = subprocessLauncher(cmdDagcon, outDagcon)
	outDagcon.close()

def daccord(suffix):
	outDaccord = open("corrected_by_daccord" + suffix + ".fa", 'w')
	cmdDaccord = "daccord db_simulatedLR" + suffix + ".las db_simulatedLR" + suffix + ".db"
	p = subprocessLauncher(cmdDaccord, outDaccord)
	outDaccord.close()


def proovread(suffix, currentDirectory):
	cmdProovread =  "proovread -l simulatedLR" + suffix + ".fa -s simulatedSR" + suffix + ".fa --pre " + currentDirectory + "/outProovread"
	p = subprocessLauncher(cmdProovread)
	cmdfa = "fastq2fasta.sh -q " + currentDirectory + "/outProovread/outProovread.untrimmed.fq -f " + currentDirectory + "/corrected_by_Proovread" + suffix + ".fa"
	p = subprocessLauncher(cmdfa)
	#cmdfa = "/home/marchet/scripts/fastq2fasta.sh -q outProovread/outProovread.trimmed.fq -f /home/marchet/hybrid_correction/ASTER/benchs/corrected_by_ProovreadTrimmed" + suffix + ".fa"
	#p = subprocessLauncher(cmdfa)
	cmdrm = "rm -r outProovread"
	subprocess.check_output(['bash','-c', cmdrm])


def hgcolor(suffix):
	cmdmkdir = "mkdir tmp"
	subprocess.check_output(['bash','-c', cmdmkdir])
	cmdHgc = "HG-CoLoR --longreads simulatedLR" + suffix + ".fa --shortreads simulatedSR" + suffix + ".faS --out corrected_by_HGColor" + suffix + ".fa --tmpdir tmp"
	p = subprocessLauncher(cmdHgc)


def launchCorrectors(currentDirectory, listCorrectors = ["LoRDEC", "colorMap", "LoRMA", "MECAT", "PBDagCon", "daccord", "Proovread", "msa_exon","msa_isoform","msa_sparc"], skipped=["10","50","100"], abund=[ "90", "75 ","50"], nbSR = 100, nbLR = 10):
	for sizeSkipped in skipped:
	#~ for sizeSkipped in ["10", "50", "100"]:
		for relAbund in abund:
	#~ for sizeSkipped in ["100"]:
		#~ for relAbund in [ "50"]:
			suffix = "_size_" + sizeSkipped + "_abund_" + relAbund
			# simulation
			#~ simulation(sizeSkipped, relAbund, suffix, nbSR, nbLR)
			################## LORDEC ##########################
			if "LoRDEC" in listCorrectors:
				lordec(suffix)
			#~ ################## COLORMAP ##########################
			if "colorMap" in listCorrectors:
				colormap(suffix)
			#~ ################## LORMA ##########################
			if "LoRMA" in listCorrectors:
				lorma(suffix)
			#~ ################## MECAT ##########################
			if "MECAT" in listCorrectors:
				mecat(suffix)
			#~ ################## HG-COLOR ##########################
			if "HGColor" in listCorrectors:
				hgcolor(sufix)
			################## DAZZDB + DALIGNER ##########################
			if "PBDagCon" in listCorrectors or "daccord" in listCorrectors:
			#~ #for daccord and pbdagcon: generate subreads alignment with daligner + db with dazz
			#~ #make pacbio headers for reads
				dazzdb_daligner(suffix)
			################## PBDAGCON ##########################
				if "PBDagCon" in listCorrectors:
				#~ #pbdagcon with subreads generated by daligner
					pbdagcon(suffix)
			################## DACCORD ##########################
			#daccord with subreads generated by daligner
				if "daccord" in listCorrectors:
					daccord(suffix)
				cmdrm = "rm *.las; rm db*"
				try:
					subprocess.check_output(['bash','-c', cmdrm])
				except subprocess.CalledProcessError:
					pass
			################## PROOVREAD ##########################
			if "Proovread" in listCorrectors:
				proovread(suffix, currentDirectory)
			################## MSA ##########################
			if "msa_exon" in listCorrectors:
				msa(suffix, "msa_exon")
			if "msa_isoform" in listCorrectors:
				msa(suffix, "msa_isoform")
			if "msa_sparc" in listCorrectors:
				msa(suffix, "msa_sparc")

	cmdrm = "rm *.h5 ; rm candidates*  ; rm trashme*"
	try:
		subprocess.check_output(['bash','-c', cmdrm])
	except subprocess.CalledProcessError:
		pass


# check results by aligning on the reference sequences
def alignOnRef(listCorrectors = ["LoRDEC", "colorMap", "LoRMA", "MECAT", "PBDagCon", "daccord", "Proovread","msa_exon","msa_isoform","msa_sparc"],skipped=[10,50,100], abund=[ 90, 75 ,50]):
	for soft in listCorrectors:
		for sizeSkipped in skipped:
		#~ for sizeSkipped in [100, 50 ,10]:
			for relAbund in abund:
		#~ for sizeSkipped in [100]:
			#~ for relAbund in [50]:
				suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
				samFile = open("results" + soft + suffix + ".sam", 'w')
				cmdAlign = "./Complete-Striped-Smith-Waterman-Library/src/ssw_test refSequences" + suffix + ".fa corrected_by_" + soft + suffix + ".fa  -c -s"
				p = subprocessLauncher(cmdAlign, samFile)






def alignOnRefMsa(skipped, abund, soft):
	switches = []
	#~ soft = "MSA"
	
	for sizeSkipped in skipped:
	#~ for sizeSkipped in [100, 50 ,10]:
		for relAbund in abund:
		#~ for relAbund in [50, 75, 90]:
			suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
			fileNameC0 = "corrected_by_" + soft + "_c0" + suffix + ".fa"
			fileNameC1 = "corrected_by_" +soft + "_c1" + suffix + ".fa"
			isoform = None
			clusterToIsoform = dict()
			if (os.path.isfile(fileNameC0) and os.access(fileNameC0, os.R_OK)) or (os.path.isfile(fileNameC1) and os.access(fileNameC1, os.R_OK)):
				c0 = open(fileNameC0, 'r')
				c0Lines = c0.readlines()
				prec = None
				switch = False
				for r in c0Lines:
					if ">" in r:
						isoform = r.split(">")[1].split("_")[0]
						if prec is not None:
							if isoform != prec:
								switch = True
						else:
							prec = isoform
				if not switch:
					clusterToIsoform[fileNameC0] = isoform
				c1 = open(fileNameC1, 'r')
				c1Lines = c1.readlines()
				prec = None
				for r in c1Lines:
					if ">" in r:
						isoform = r.split(">")[1].split("_")[0]
						if prec is not None:
							if isoform != prec:
								switch = True
						else:
							prec = isoform
				if switch:
					alignOnRef([soft])
					switch = False
					#~ return 0
					
				else:
					clusterToIsoform[fileNameC1] = isoform
					for fileNameClust in clusterToIsoform.keys():
						isoform = clusterToIsoform[fileNameClust]
						cmdGrep = "grep "+ isoform + " refSequences" + suffix + ".fa -A 1 > refSequences" + isoform + suffix + ".fa"
						subprocess.check_output(['bash','-c', cmdGrep])
						cmdGrep = "grep "+ isoform + " " + fileNameClust + " -A 1 > toalign.fa"
						subprocess.check_output(['bash','-c', cmdGrep])
						samFile = open("results" + isoform + soft + suffix + ".sam", 'w')
						#~ cmdAlign = "./Complete-Striped-Smith-Waterman-Library/src/ssw_test refSequences" + isoform + suffix + ".fa " + fileNameClust + " -c -s"
						cmdAlign = "./Complete-Striped-Smith-Waterman-Library/src/ssw_test refSequences" + isoform + suffix + ".fa toalign.fa -c -s"
						p = subprocessLauncher(cmdAlign, samFile)
						samFile.close()
					
					
					
				cmdCat = "cat resultsinclusion" + soft + suffix + ".sam resultsexclusion" + soft + suffix + ".sam > results" + soft + suffix + ".sam"
				subprocess.check_output(['bash','-c', cmdCat])
			else:
				#~ alignOnRef(["MSA"])
				alignOnRef([soft])
				#~ return 0
				switch = False
			switches.append(switch)
	return switches





def getExpectedLength(suffix):
	expectedLengths = {"exclusion":0, "inclusion":0}
	refFile = open("refSequences" + suffix + ".fa", 'r')
	lines = refFile.readlines()
	for l in lines:
		if ">" in l:
			targetType = l[1:-1]
		else:
			expectedLengths[targetType] = len(l) - 1
	return expectedLengths



def readSam(soft, suffix):
	blockResults = dict()
	lenResults = dict()
	pathSam = "results" + soft + suffix + ".sam"
	
	if os.path.exists(pathSam) and os.path.getsize(pathSam) > 0:
		#~ print("%%%%%%%%%%%%%%%%%%%%%%%ù", pathSam)
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
						#~ print(b,r)
						gapsLength += int(r)
			if query not in lenResults:
				lenResults[query] = {target:[length, alnLength -start, gapsLength]}
			else:
				lenResults[query][target] = [length, alnLength -start, gapsLength]
		#~ print(query, target)
		return ( start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries)
	else:
		return (None,) * 8




def getIsoform(blockResults, lenResults, suffix, queries, nbIncl, nbExcl, soft):
	countIncl = 0
	countExcl = 0
	countAssigned = 0
	countInverted = 0
	outI = open("corrected_reads_inclusion" + suffix + ".fa", 'w')
	outE = open("corrected_reads_exclusion" + suffix + ".fa", 'w')
	meanSizes = {"exclusion":{"realSize":[], "alignedSize":[]}, "inclusion":{"realSize":[], "alignedSize":[]}}
	for querySeq, aln in blockResults.items():
		#~ print("1", querySeq, aln)
		if len(aln) == 1 or len(lenResults[querySeq]) == 1: # corrected read could be aligned in 1 block on one ref transcript
			nill, targetType = next(iter(aln.items()))
		
		else: #in case sequences are not aligned in one block, keep the reference linked with the smaller cumulative length of gaps
			minGapsLen = None
			targetType = None
			for tar in lenResults[querySeq]:
				if minGapsLen is None:
					minGapsLen = lenResults[querySeq][tar][-1]
					targetType = tar
				else:
					if minGapsLen > lenResults[querySeq][tar][-1]:
						minGapsLen = lenResults[querySeq][tar][-1]
						targetType = tar
		#~ print("2", targetType)
		toWRatio = str(nbIncl)
		if  targetType == "inclusion":
			#~ print("target type inclusion")
			countAssigned += 1
			countIncl += 1
			meanSizes[targetType]["realSize"].append(lenResults[querySeq][targetType][0])
			meanSizes[targetType]["alignedSize"].append(lenResults[querySeq][targetType][1])
			#~ ratio = round(meanSizes[targetType]["realSize"][-1]/expectedLengths["inclusion"],3) *100
			#~ out2.write(soft + " inclusion " +  str(relAbund)  + " " + str(ratio) + " " + str(sizeSkipped) +"\n")
			#~ if soft == "MSA":
			if "msa" in soft:
				outI.write(">" +  querySeq + "\n" + queries[querySeq] + "\n")
			else:
				outI.write(">" +  querySeq.split('_')[1] + "\n" + queries[querySeq] + "\n")
			if querySeq.split('_')[0] != "inclusion":
				#~ print("happens", targetType, querySeq)
				countInverted += 1
		toWRatio += " " + str(countInverted)

		count = 0
		if  targetType == "exclusion":
			#~ print("target type exclusion")
			countAssigned += 1
			countExcl += 1
			meanSizes[targetType]["realSize"].append(lenResults[querySeq][targetType][0])
			meanSizes[targetType]["alignedSize"].append(lenResults[querySeq][targetType][1])
			#~ ratio = round(meanSizes[targetType]["realSize"][-1]/expectedLengths["exclusion"],3)*100 
			#~ out2.write(soft + " exclusion " +  str(relAbund)  + " " + str(ratio) + " " + str(sizeSkipped) +"\n")
			if "msa" in soft:
			#~ if soft == "MSA":
				outE.write(">" +  querySeq + "\n" + queries[querySeq] + "\n")
			else:
				outE.write(">" +  querySeq.split('_')[1] + "\n" + queries[querySeq] + "\n")
			if querySeq.split('_')[0] != "exclusion":
				#~ print("happens", targetType, querySeq)
				countInverted += 1
				count += 1
		toWRatio +=  " " + str(count)  + " " + str(nbExcl)
		
	outRatio = open(soft + "_ratio_" + suffix + ".txt", 'w')
	#nb_incl nb_incl_corrected_excl nb_excl nb_excl_corrected_incl
	outRatio.write(toWRatio)
	outRatio.close()
	return meanSizes, countIncl, countExcl, countAssigned, countInverted





def getFileReadNumber(fileName):
	cmdGrep = """grep ">" -c """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return int(val.decode('ascii'))



def computeResults(listCorrectors = ["LoRDEC", "colorMap", "LoRMA", "MECAT", "PBDagCon", "daccord", "Proovread","msa_exon","msa_isoform","msa_sparc"], skipped=[10,50,100], abund=[ 90, 75 ,50]):
	out = open("corrected_to_inclusion.txt", 'w')
	
	out2 = open("corrected_sizes.txt", 'w')
	out.write("corrector size_skipped_exon ratio percent\n")
	out2.write("corrector isoform ratio percent size_skipped_exon\n")
	outSize = open("sizes_reads.txt", 'w')
	outSize.write("soft size skipped abund\n")
	for soft in listCorrectors:
	
			
		for sizeSkipped in skipped:
		#~ for sizeSkipped in [10,50,100]:
			for relAbund in abund:
			#~ for relAbund in [ 90, 75 ,50]:
		#~ for sizeSkipped in [100]:
			#~ for relAbund in [50]:
				
				suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)

				expectedLengths = getExpectedLength(suffix)
				
				start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries = readSam(soft, suffix)
				if start is not None:
					nbIncl = getFileReadNumber("perfect_reads_inclusion" + suffix + ".fa")
					nbExcl = getFileReadNumber("perfect_reads_exclusion" + suffix + ".fa")
					meanSizes, countIncl, countExcl, countAssigned, countInverted = getIsoform(blockResults, lenResults, suffix, queries, nbIncl, nbExcl, soft)
					readNumber = len(blockResults.keys())
					meanReadsSize = round(sum(readsSize)*1.0/len(readsSize),2) if len(readsSize) > 0 else 0
					meanInclusionCorrectedSize = round(sum(meanSizes["inclusion"]["alignedSize"])*1.0/len(meanSizes["inclusion"]["alignedSize"]),2) if len(meanSizes["inclusion"]["alignedSize"]) > 0 else 0
					meanExclusionCorrectedSize = round(sum(meanSizes["exclusion"]["alignedSize"])*1.0/len(meanSizes["exclusion"]["alignedSize"]),2) if len(meanSizes["exclusion"]["alignedSize"]) > 0 else 0

					ratioLenI = round(meanInclusionCorrectedSize*100/expectedLengths["inclusion"],2)
					ratioLenE = round(meanExclusionCorrectedSize*100/expectedLengths["exclusion"],2)
					
					print("#############################")
					print(soft, "correction for inclusion size", sizeSkipped, ", ratio:", relAbund, "/", 100-relAbund)
					print("Corrected inclusion reads in output:", round(countIncl*100.0/(countIncl+countExcl),2) if countIncl+countExcl != 0 else 0, "%")
					print("Corrected exclusion reads in output:", round(countExcl*100.0/(countExcl+countIncl),2) if countIncl+countExcl != 0 else 0, "%")
					print("Mean length corrected reads:", meanReadsSize, ", mean length of inclusion reads:", meanInclusionCorrectedSize, "(", ratioLenI, "% of real size ), mean length of exclusion reads:", meanExclusionCorrectedSize, "(", ratioLenE, "% of real size )")
					outSize.write(soft + " " + str(ratioLenI) + " " +  str(sizeSkipped) + " "+ str(relAbund) + "\n")
					outSize.write(soft + " " + str(ratioLenE) + " " +  str(sizeSkipped) + " "+ str(relAbund) +"\n")
					correctedToIncl = round(countIncl*100.0/(countIncl+countExcl),2) if countIncl+countExcl != 0 else 0
					if not correctedToIncl == 0 :
						out.write(soft + " " + str(sizeSkipped) + " " + str(relAbund) + " " + str(correctedToIncl) +"\n")
						out2.write(soft + " inclusion " +  str(relAbund)  + " " + str(ratioLenI) + " " + str(sizeSkipped) +"\n")
						out2.write(soft + " exclusion " +  str(relAbund)  + " " + str(ratioLenE) + " " + str(sizeSkipped) +"\n")

					#~ ### launch corrector benchmark"
					if "msa" not in soft :
					#~ if soft != "MSA":
						if countExcl > 0:
							cmdBench = "python3 ./benchmark-long-read-correction/benchmark.py -c corrected_reads_exclusion" + suffix + ".fa -u uncorrected_reads_exclusion" + suffix + ".fa -r perfect_reads_exclusion" + suffix + ".fa"
							print("Exclusion")
							subprocessLauncher(cmdBench)
							cmdMv = "mv outPositions.txt outPositionsExclusion_" + soft + suffix + ".fa"
							subprocess.check_output(['bash','-c', cmdMv])

							cmdSed = "sed -i 's/unknown/" + soft + "/g' precision.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " precision.txt >> precision_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
							
							cmdSed = "sed -i 's/unknown/" + soft + "/g' recall.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " recall.txt >> recall_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
							
							cmdSed = "sed -i 's/unknown/" + soft + "/g' correct_base_rate.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " correct_base_rate.txt >> correct_base_rate_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
							
						if countIncl > 0:
							cmdBench = "python3 ./benchmark-long-read-correction/benchmark.py -c corrected_reads_inclusion" + suffix + ".fa -u uncorrected_reads_inclusion" + suffix + ".fa -r perfect_reads_inclusion" + suffix + ".fa"
							print("Inclusion")
							subprocessLauncher(cmdBench)
							cmdMv = "mv outPositions.txt outPositionsInclusion_" + soft + suffix + ".fa"
							subprocess.check_output(['bash','-c', cmdMv])

							cmdSed = "sed -i 's/unknown/" + soft + "/g' precision.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " precision.txt >> precision_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
							
							cmdSed = "sed -i 's/unknown/" + soft + "/g' recall.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " recall.txt >> recall_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])

							cmdSed = "sed -i 's/unknown/" + soft + "/g' correct_base_rate.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " correct_base_rate.txt >> correct_base_rate_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
					else:
						print("Exclusion")
						if countExcl > 0:
							cmdHead = "head -2 corrected_reads_exclusion" + suffix + ".fa > corrected.fa"
							subprocess.check_output(['bash','-c', cmdHead])
							cmdHead = "head -2 uncorrected_reads_exclusion" + suffix + ".fa > uncorrected.fa"
							subprocess.check_output(['bash','-c', cmdHead])
							cmdHead = "head -2 perfect_reads_exclusion" + suffix + ".fa > perfect.fa"
							subprocess.check_output(['bash','-c', cmdHead])
							
							cmdBench = "python3 ./benchmark-long-read-correction/benchmark.py -c corrected.fa -u uncorrected.fa -r perfect.fa"
							
							subprocessLauncher(cmdBench)
							cmdMv = "mv outPositions.txt outPositionsExclusion_" + soft +  suffix + ".fa"
							subprocess.check_output(['bash','-c', cmdMv])

							cmdSed = "sed -i 's/unknown/" + soft + "/g' precision.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " precision.txt >> precision_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
							
							cmdSed = "sed -i 's/unknown/" + soft + "/g' recall.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " recall.txt >> recall_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])

							cmdSed = "sed -i 's/unknown/" + soft + "/g' correct_base_rate.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " correct_base_rate.txt >> correct_base_rate_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
						else:
							print("Recall 0 Precision 0 correct_bases_rate 0")

						print("Inclusion")
						if countIncl > 0:
							
							cmdHead = "head -2 corrected_reads_inclusion" + suffix + ".fa > corrected.fa"
							subprocess.check_output(['bash','-c', cmdHead])
							cmdHead = "head -2 uncorrected_reads_inclusion" + suffix + ".fa > uncorrected.fa"
							subprocess.check_output(['bash','-c', cmdHead])
							cmdHead = "head -2 perfect_reads_inclusion" + suffix + ".fa > perfect.fa"
							subprocess.check_output(['bash','-c', cmdHead])
							
							cmdBench = "python3 ./benchmark-long-read-correction/benchmark.py -c corrected.fa -u uncorrected.fa -r perfect.fa"
							
							subprocessLauncher(cmdBench)
							cmdMv = "mv outPositions.txt outPositionsInclusion_" + soft + suffix + ".fa"
							subprocess.check_output(['bash','-c', cmdMv])
							
							cmdSed = "sed -i 's/unknown/" + soft + "/g' precision.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " precision.txt >> precision_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
							
							cmdSed = "sed -i 's/unknown/" + soft + "/g' recall.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " recall.txt >> recall_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])

							cmdSed = "sed -i 's/unknown/" + soft + "/g' correct_base_rate.txt"
							subprocess.check_output(['bash','-c', cmdSed])
							cmdCat = "grep " + soft + " correct_base_rate.txt >> correct_base_rate_tmp.txt"
							subprocess.check_output(['bash','-c', cmdCat])
						else:
							print("Recall 0 Precision 0 correct_bases_rate 0")
	cmdMv = "mv recall_tmp.txt recall.txt"
	subprocess.check_output(['bash','-c', cmdMv])
	cmdMv = "mv precision_tmp.txt precision.txt"
	subprocess.check_output(['bash','-c', cmdMv])

	cmdMv = "mv correct_base_rate_tmp.txt correct_base_rate.txt"
	subprocess.check_output(['bash','-c', cmdMv])

#todo mean between inclusion and exclusion in R



def simulateReads(covSR, covLR, skipped, abund):
	for sizeSkipped in skipped:
		for relAbund in  abund:
	#~ for sizeSkipped in [ "100"]:
		#~ for relAbund in ["50"]:
			suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
			# simulation
			simulation(str(sizeSkipped), str(relAbund), suffix, covSR, covLR)

def plotResults(skipped, abund):
	cmdR = "Rscript corrected_to_major.R corrected_to_inclusion.txt"
	subprocessLauncher(cmdR)
	cmdR = "Rscript corrected_to_major_self.R corrected_to_inclusion.txt"
	subprocessLauncher(cmdR)
	
	cmdR = "Rscript plot_recall.R recall.txt"
	subprocessLauncher(cmdR)
	cmdR = "Rscript plot_precision.R precision.txt"
	subprocessLauncher(cmdR)

	cmdR = "Rscript plot_correct_base_rate.R correct_base_rate.txt"
	subprocessLauncher(cmdR)

	cmdR =  "Rscript plot_size.R sizes_reads.txt "
	subprocessLauncher(cmdR)


	for sizeSkipped in skipped:
	#~ for sizeSkipped in [10,50,100]:
			for relAbund in abund:
		#~ for sizeSkipped in [100]:
			#~ for relAbund in [50]:
				
				suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund)
				cmdR = "Rscript error_frequencies_inclusion.R outPositionsInclusion_MSA" + suffix + ".fa 800 " + str(800+sizeSkipped)
				subprocessLauncher(cmdR)
				try:
					cmdMv = "mv error_frequencies_by_positions.png error_frequencies_by_positionsInclusion"  + suffix + ".png"
					subprocess.check_output(['bash','-c', cmdMv])
				except subprocess.CalledProcessError:
					pass
				cmdR = "Rscript error_frequencies_exclusion.R outPositionsExclusion_MSA" + suffix + ".fa 800"
				subprocessLauncher(cmdR)
				try:
					cmdMv = "mv error_frequencies_by_positions.png error_frequencies_by_positionsExclusion"  + suffix + ".png"
					subprocess.check_output(['bash','-c', cmdMv])
				except subprocess.CalledProcessError:
					pass



def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))
	#~ listCorrectors = ["LoRDEC", "colorMap", "LoRMA", "MECAT", "PBDagCon", "daccord", "Proovread"]
	#~ listCorrectors = ["MSA","daccord"]
	#~ listCorrectors = ["MSA","daccord",  "LoRMA", "MECAT", "PBDagCon"]
	#~ listCorrectors = ["LoRDEC", "colorMap", "LoRMA", "MECAT", "PBDagCon", "daccord"]
	#~ listCorrectors = ["msa_exon", "msa_isoform", "msa_sparc"]
	#~ listCorrectors = ["msa_exon"]
	#~ listCorrectors = ["msa_isoform"]
	listCorrectors = ["msa_exon", "msa_isoform"]
	#~ listCorrectors = ["Proovread"]
	covSR = 1
	covLR = 20
	#~ skipped = [100]
	#~ abund = [50]
	skipped = [100,50]
	abund = [75,50,90]
	skippedS = [str(r) for r in skipped]
	abundS = [str(r) for r in abund]
	readsFileName = None
	if len(sys.argv) > 1:
		correctors = sys.argv[1]
		if correctors != "all":
			listCorrectors = [correctors]
		if len(sys.argv) > 3:
			readsFileName = argv[3]
			if len(sys.argv) > 4:
				covSR = int(sys.argv[4])
				covLR = int(sys.argv[5])
	if readsFileName is None:
		simulateReads(covSR, covLR, skipped, abund)
	launchCorrectors(currentDirectory,listCorrectors, skippedS, abundS, covSR, covLR)

	listCorrectorsCpy = copy.copy(listCorrectors)
	msaIn = []
	for msa in ["msa_exon", "msa_isoform", "msa_sparc"]:
		
		if  msa in listCorrectors:
			listCorrectorsCpy.remove(msa)
			msaIn.append(msa)
	if len(msaIn) !=0:
		for soft in msaIn:
			alignOnRef(listCorrectorsCpy, skipped, abund)
			alignOnRefMsa(skipped, abund, soft)
	else:
		alignOnRef(listCorrectors, skipped, abund)
	computeResults(listCorrectors, skipped, abund)
	plotResults(skipped, abund)

if __name__ == '__main__':
	main()
