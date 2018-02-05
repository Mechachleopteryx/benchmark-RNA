
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
import benchmark_msa
import benchmark_self_correctors
import utils


		
def lordec(suffix):
	cmdLordec = "lordec-correct -i simulatedLR" + suffix + ".fa -2 simulatedSR" + suffix + ".fa -o corrected_by_LoRDEC" + suffix + ".fa -k 21 -s 2 -T 4"
	p = subprocessLauncher(cmdLordec)


def colormap(suffix):
	cmdColorMap = "runCorr.sh simulatedLR"  + suffix + ".fa simulatedSR"  + suffix + ".fa colorMapDir pre 4"
	p = subprocessLauncher(cmdColorMap)
	cmdmv = "mv colorMapDir/pre_sp.fasta corrected_by_colorMap" + suffix + ".fa"
	subprocess.check_output(['bash','-c', cmdmv])


def lorma(suffix, currentDirectory):
	cmdcp = "cp " + currentDirectory + "/simulatedLR" + suffix + ".fa " + currentDirectory + "/copy.fa"
	subprocess.check_output(['bash','-c', cmdcp])
	cmdLorma = "/home/marchet/bin/LoRMA-0.4/build/lorma.sh -threads 4 " + currentDirectory + "/copy.fa -s -n"
	p = subprocessLauncher(cmdLorma)
	try:
		if checkIfFile(currentDirectory + "/trim.fasta"):
			cmdmv = "/home/marchet/scripts/convertOneLineFasta.py " + currentDirectory + "/trim.fasta > " + currentDirectory + "/corrected_by_LoRMA" + suffix + ".fa"
			subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	try:
		if checkIfFile(currentDirectory + "/discarded.fasta"):
			cmdmv = "/home/marchet/scripts/convertOneLineFasta.py " + currentDirectory + "/discarded.fasta > " + currentDirectory + "/discarded_by_LoRMA" + suffix + ".fa"
			subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	try:
		if checkIfFile(currentDirectory + "/final.fasta"):
			cmdmv = "/home/marchet/scripts/convertOneLineFasta.py " + currentDirectory + "/final.fasta > " + currentDirectory + "/corrected_by_LoRMA" + suffix + ".fa"
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





def launchCorrector(currentDirectory, corrector , sizeSkipped, relAbund, coverage, error):
	#~ for sizeSkipped in skipped:
	#~ for sizeSkipped in ["10", "50", "100"]:
		#~ for relAbund in abund:
	#~ for sizeSkipped in ["100"]:
		#~ for relAbund in [ "50"]:
	suffix = "_size_" + sizeSkipped + "_abund_" + relAbund + "_cov_" + str(coverage) + "_err_" + str(error)
	# simulation
	#~ simulation(sizeSkipped, relAbund, suffix, nbSR, nbLR)
	################## LORDEC ##########################
	if "LoRDEC" == corrector:
		lordec(suffix)
	#~ ################## COLORMAP ##########################
	if "colorMap" == corrector:
		colormap(suffix)
	#~ ################## LORMA ##########################
	if "LoRMA" == corrector:
		lorma(suffix, currentDirectory)
	#~ ################## MECAT ##########################
	if "MECAT" == corrector:
		mecat(suffix)
	#~ ################## HG-COLOR ##########################
	if "HGColor" == corrector:
		hgcolor(sufix)
	################## DAZZDB + DALIGNER ##########################
	if "PBDagCon" == corrector or "daccord" == corrector:
	#~ #for daccord and pbdagcon: generate subreads alignment with daligner + db with dazz
	#~ #make pacbio headers for reads
		dazzdb_daligner(suffix)
	################## PBDAGCON ##########################
		if "PBDagCon" == corrector:
		#~ #pbdagcon with subreads generated by daligner
			pbdagcon(suffix)
	################## DACCORD ##########################
	#daccord with subreads generated by daligner
		if "daccord" == corrector:
			daccord(suffix)
		cmdrm = "rm *.las; rm db*"
		try:
			subprocess.check_output(['bash','-c', cmdrm])
		except subprocess.CalledProcessError:
			pass
	################## PROOVREAD ##########################
	if "Proovread" == corrector:
		proovread(suffix, currentDirectory)

	#~ cmdrm = "rm *.h5 ; rm candidates*  ; rm trashme*"
	#~ try:
		#~ subprocess.check_output(['bash','-c', cmdrm])
	#~ except subprocess.CalledProcessError:
		#~ pass


	


# check results by aligning on the reference sequences
def alignOnRef(soft,sizeSkipped, relAbund, currentDirectory, coverage, error):
	suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund) + "_cov_" + str(coverage) + "_err_" + str(error)
	#~ listFileNames = getFiles(resultDirectory, "corrected_by_" + soft + "*.fa")
	#~ for fileC in listFileNames:
	samFile = open("results" + soft + suffix + ".sam", 'w')
	cmdAlign = "/home/marchet/bin/Complete-Striped-Smith-Waterman-Library/src/ssw_test " + currentDirectory + "/refSequences" + suffix + ".fa "+ currentDirectory + "/corrected_by_" + soft + suffix + ".fa  -c -s"
	p = subprocessLauncher(cmdAlign, samFile)

def makeCorrectedHeadersAndCompareFromSam(refIsoformTypesToCounts, blockResults):

	for querySeq, aln in blockResults.items():
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
		toWRatio = str(nbIncl)
		if  targetType == "inclusion":
			countAssigned += 1
			countIncl += 1
			meanSizes[targetType]["realSize"].append(lenResults[querySeq][targetType][0])
			meanSizes[targetType]["alignedSize"].append(lenResults[querySeq][targetType][1])
			outI.write(">" +  querySeq.split('_')[1] + "\n" + queries[querySeq] + "\n")
			if querySeq.split('_')[0] != "inclusion":
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



def computeResults(soft,sizeSkipped, relAbund, currentDirectory, coverage, error, refIsoformTypesToSeq):					
	suffix = "_size_" + str(sizeSkipped) + "_abund_" + str(relAbund) + "_cov_" + str(coverage) + "_err_" + str(error)

	
	#~ for isofType in refIsoformTypesToSeq:
		#~ expectedLengths = getPerfectSequenceLength(currentDirectory + "/perfect_reads_" + isofType + suffix + ".fa")
	pathSam  = currentDirectory + "/results" + soft + suffix + ".sam"
	start, readsSize, resultAln, gapsLength, blockResults, alnLength, lenResults, queries = readSam(soft, suffix, pathSam, currentDirectory)
	if start is not None:
		nbIncl = getFileReadNumber("perfect_reads_inclusion" + suffix + ".fa")
		nbExcl = getFileReadNumber("perfect_reads_exclusion" + suffix + ".fa")
		meanSizes, countIncl, countExcl, countAssigned, countInverted = getIsoform(blockResults, lenResults, suffix, queries, nbIncl, nbExcl, soft)
		readNumber = len(blockResults.keys())
		meanReadsSize = round(sum(readsSize)*1.0/len(readsSize),2) if len(readsSize) > 0 else 0
		meanInclusionCorrectedSize = round(sum(meanSizes["inclusion"]["alignedSize"])*1.0/len(meanSizes["inclusion"]["alignedSize"]),2) if len(meanSizes["inclusion"]["alignedSize"]) > 0 else 0
		meanExclusionCorrectedSize = round(sum(meanSizes["exclusion"]["alignedSize"])*1.0/len(meanSizes["exclusion"]["alignedSize"]),2) if len(meanSizes["exclusion"]["alignedSize"]) > 0 else 0

	#~ ratioLenI = round(meanInclusionCorrectedSize*100/expectedLengths["inclusion"],2)
	#~ ratioLenE = round(meanExclusionCorrectedSize*100/expectedLengths["exclusion"],2)

	#~ print("#############################")
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
