#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>


using namespace std;



char randomNucleotide(){
	switch (rand() % 4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}


char replaceRandomNucleotide(char c){
	switch (rand()%4){
		case 0:
				if(c!='A'){
						return 'A';
				}
				return replaceRandomNucleotide(c);
		case 1:
				if(c!='C'){
						return 'C';
				}
				return replaceRandomNucleotide(c);
		case 2:
				if(c!='G'){
						return 'G';
				}
				return replaceRandomNucleotide(c);
		case 3:
				if(c!='T'){
						return 'T';
				}
				return replaceRandomNucleotide(c);
	}
	return replaceRandomNucleotide(c);
}




string randomSequence(const uint length){
	string result(length, 'A');
	for(uint i(0); i < length; ++i){
		result[i] = randomNucleotide();
	}	
	return result;
}


void insertion(double rate, string& result){
	uint dice(rand() % 100);
	if(dice < rate){
		char newNucleotide(randomNucleotide());
		result.push_back(newNucleotide);
		insertion(rate, result);
	}
}


string mutateSequence(const string& referenceSequence, uint maxMutRate=13, vector <double> ratioMutation={0.37,0.09,0.54}){
	string result;
	result.reserve(5 * referenceSequence.size());
	for(uint i(0); i < referenceSequence.size(); ++i){
		uint mutRate(maxMutRate);
		double substitutionRate(mutRate * ratioMutation[0]);
		double insertionRate(mutRate * ratioMutation[1]);
		double deletionRate(mutRate * ratioMutation[2]);
		uint dice(rand() % 100);
		if (dice <substitutionRate ){
			//SUBSTITUTION
			char newNucleotide(replaceRandomNucleotide(referenceSequence[i]));
			result.push_back(newNucleotide);
			continue;
		} else if(dice < deletionRate+substitutionRate){
			//DELETION
			uint dice2(rand() % 100);
			while (dice2 < deletionRate+substitutionRate){ // deletions larger than 1
				++i;
				dice2 = rand() % 100;
			}
			continue;
		} else if (dice < deletionRate + substitutionRate + insertionRate){
			//INSERTION
			char newNucleotide(randomNucleotide());
			result.push_back(referenceSequence[i]);
			result.push_back(newNucleotide);
			insertion(deletionRate + substitutionRate + insertionRate, result); // larger than 1 insertions
			continue;
		} else {
			result.push_back(referenceSequence[i]);
		}
	}
	return result;
}



vector<string> generateAlternativeTranscriptReferences(uint sizeSkippedExon, uint sizeExons=200, uint transcriptNumber=4, uint exonNumber=8){
	vector<string> result;
	vector<string> exonList;
	// generate sequences for regular exons
	for(uint i(0); i < exonNumber; ++i){
		exonList.push_back(randomSequence(sizeExons));
	}
	// skipped exons
	exonList.push_back(randomSequence(sizeSkippedExon));
	exonList.push_back(randomSequence(sizeSkippedExon)); // 2 alt exons
	// transcripts
	for(uint i(0); i < transcriptNumber; ++i){
		string transcript;
		transcript.reserve((exonNumber+2)*sizeExons);
		if (i < 3){  // inclusion  isoform
			uint ii(0);
			while (ii < exonList.size()/2){
				transcript += exonList[ii];
				++ii;
			}
			if (i==0){ //inclusion x2  isoform
				transcript += exonList.back(); // included exon in the middle
				transcript += exonList[exonList.size() - 2]; // included exon in the middle
			} else if (i==1){
				transcript += exonList.back(); // included exon in the middle
			} else {
				transcript += exonList[exonList.size() - 2]; 
			}
			while (ii < exonList.size() - 2){
				transcript += exonList[ii];
				++ii;
			}
		}else{ // exclusion isoform
			for(uint ii(0); ii < exonNumber; ++ii){
				transcript += exonList[ii];
			}
		}
		// polyA tails
		for (uint a(0); a < 50; ++a){
			transcript += 'A';
		}
		result.push_back(transcript);
	}
	return result; // [0]inclusion, [1]exclusion
}






void generateReads(uint inclusionAbundance, vector<string>& transcripts, string& outSuffix, uint lrCoverage=100, uint errorRate=13, uint srLength=150){
	ofstream outLR("simulatedLR" + outSuffix + ".fa");
	//~ ofstream outSR("simulatedSR" + outSuffix + ".fa");
	ofstream outRef("refSequences" + outSuffix + ".fa");
	ofstream outPerfectIncl12("perfect_reads_inclusion12" + outSuffix + ".fa");
	ofstream outPerfectIncl1("perfect_reads_inclusion1" + outSuffix + ".fa");
	ofstream outPerfectIncl2("perfect_reads_inclusion2" + outSuffix + ".fa");
	ofstream outPerfectExcl("perfect_reads_exclusion" + outSuffix + ".fa");
	ofstream outUncoIncl12("uncorrected_reads_inclusion12" + outSuffix + ".fa");
	ofstream outUncoIncl1("uncorrected_reads_inclusion1" + outSuffix + ".fa");
	ofstream outUncoIncl2("uncorrected_reads_inclusion2" + outSuffix + ".fa");
	ofstream outUncoExcl("uncorrected_reads_exclusion" + outSuffix + ".fa");
	string refSequence;
	uint numberLRInclusion12(inclusionAbundance * lrCoverage / 100 );
	uint numberLRInclusion1((100-inclusionAbundance) * (double)lrCoverage / (3.0*100) );
	uint numberLRInclusion2((100-inclusionAbundance) * (double)lrCoverage / (3.0*100) );
	uint numberLRExclusion((100-inclusionAbundance) * (double)lrCoverage / (3.0*100) );
	//~ cout << (double)lrCoverage / (3.0*100) << endl;
	//~ cout << "number of reads " <<  numberLRInclusion12 << " " << numberLRInclusion1 << " " << numberLRInclusion2 << " " << numberLRExclusion << endl;
	//~ uint numberSRInclusion(5);
	//~ uint numberSRExclusion(6);
	//~ uint numberLRInclusion(900);
	//~ uint numberLRExclusion((1000);
	//~ cout << numberLRInclusion<< " " << numberLRExclusion << endl;
	//~ cout << numberSRInclusion<< " " << numberSRExclusion << endl;
	string read;
	// LONG READS
	for (uint i(0); i < numberLRInclusion12; ++i){
		read = mutateSequence(transcripts[0], errorRate);
		outLR << ">inclusion12_" << i <<  endl << read << endl;
		outPerfectIncl12 << ">" << i << endl <<  transcripts[0] << endl;
		outUncoIncl12 << ">" << i << endl << read << endl;
	}
	for (uint i(0); i < numberLRInclusion1; ++i){
		read = mutateSequence(transcripts[1], errorRate);
		outLR << ">inclusion1_" << i <<  endl << read << endl;
		outPerfectIncl1 << ">" << i << endl <<  transcripts[1] << endl;
		outUncoIncl1 << ">" << i << endl << read << endl;
	}
	for (uint i(0); i < numberLRInclusion2; ++i){
		read = mutateSequence(transcripts[2], errorRate);
		outLR << ">inclusion2_" << i <<  endl << read << endl;
		outPerfectIncl2 << ">" << i << endl <<  transcripts[2] << endl;
		outUncoIncl2 << ">" << i << endl << read << endl;
	}
	for (uint i(0); i < numberLRExclusion; ++i){
		read = mutateSequence(transcripts[3], errorRate);
		outLR << ">exclusion_" << i << endl <<  read << endl;
		outPerfectExcl << ">" << i << endl << transcripts[3] << endl;
		outUncoExcl << ">" << i << endl << read << endl;
	}

	//~ // SHORT READS
	//~ uint start;
	//~ for (uint i(0); i < numberSRInclusion; ++i){
		//~ start = rand() % (transcripts[0].size() - srLength);
		//~ outSR << ">inclusion_" << i << endl << transcripts[0].substr(start, srLength) << endl;
	//~ }
	//~ for (uint i(0); i < numberSRExclusion; ++i){
		//~ start = rand() % (transcripts[1].size() - srLength);
		//~ outSR << ">exclusion_" << i << endl << transcripts[1].substr(start, srLength) << endl;
	//~ }

	// REF FILE
	outRef << ">inclusion12" << endl << transcripts[0] << endl;
	outRef << ">inclusion1" << endl << transcripts[1] << endl;
	outRef << ">inclusion2" << endl << transcripts[2] << endl;
	outRef << ">exclusion" << endl << transcripts[3] << endl;
}






int main(int argc, char ** argv){
	if (argc > 2){
		uint sizeSkippedExon(stoi(argv[1]));
		uint inclusionAbundance(stoi(argv[2]));
		string outSuffix("");
		uint  lrCoverage(10);
		if (argc > 3){
			outSuffix = argv[3];
			if (argc > 4){
				lrCoverage = stoi(argv[4]);
			}
		}
		srand (time(NULL));
		auto startChrono = chrono::system_clock::now();
		// generate reference sequences for alternative transcripts
		vector<string> transcripts(generateAlternativeTranscriptReferences(sizeSkippedExon));
		// generate reads (long and short)
		generateReads(inclusionAbundance, transcripts, outSuffix, lrCoverage);
		auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
		cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	} else {
		cout << "usage: ./MES_simulation <size skipped exon> <relative abudance inclusion isoform (%)> " << endl;
	}
	return 0;
}
