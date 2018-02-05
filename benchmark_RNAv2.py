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
from benchmark_msa import *
from benchmark_self_correctors import *
from utils import *












def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))

	print(currentDirectory)
	
	covSR = 1
	#~ skipped = [40,50,100]
	skipped = [100]
	#~ abund = [50,75,90,10]
	#~ abund = [50,75,90]
	#~ errorRate = [13]
	#~ errorRate = [5,9,13,15]
	errorRate = [13]
	#~ errorRate = [9,13]
	if len(errorRate) == 1:
		errorRateToKeep = errorRate[0]
	else:
		errorRateToKeep = 13

	#~ abund = [50,75,90]
	abund = [50]
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
	#~ correctors = ["msa_exon"]
	#~ correctors = ["msa_isoform"]
	#~ correctors = ["msa_exon", "msa_isoform", "msa_both"]
	correctors = ["msa_exon", "msa_both"]
	#~ correctors = ["LoRMA"]
	#~ correctors = ["msa_sparc"]
	if args.coverage is not None:
		coverage = [args.coverage]
	else:
		coverage = [20]
		#~ coverage = [10,15]
		#~ coverage = [5, 10, 15, 20, 30, 50, 100]
		#~ coverage = [10, 15, 20, 30, 50]
		#~ coverage = [10, 15, 20]
	if len(coverage) == 1:
		covForFigs = coverage[0]
	else:
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
	
	cmdFile = '''echo "soft value coverage error" > ''' + currentDirectory +  '''/all_precisions_softs.txt'''
	subprocess.check_output(['bash','-c', cmdFile])
	cmdFile = '''echo "soft value coverage error" > ''' + currentDirectory +  '''/all_recalls_softs.txt'''
	subprocess.check_output(['bash','-c', cmdFile])
	cmdFile = '''echo "soft value coverage error" > ''' + currentDirectory +  '''/all_correctRates_softs.txt'''
	subprocess.check_output(['bash','-c', cmdFile])


	simulateReads(skipped, abund, coverage, EStype, currentDirectory, errorRate)
	doBenchMsa = []
	doOthers = []
	for c in correctors:
		if "msa" in c:
			doBenchMsa.append(c)
		else:
			doOthers.append(c)
	if len(doBenchMsa) > 0:
		benchMsa(errorRate, coverage, doBenchMsa, skippedS, abundS, currentDirectory, covForFigs, errorRateToKeep, outputDirPath, outputPDFName)
	#~ if len(doOthers) > 0:
		#~ benchC
	


if __name__ == '__main__':
	main()
