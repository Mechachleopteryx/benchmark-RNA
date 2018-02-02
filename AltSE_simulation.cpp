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



vector<string> generateAlternativeTranscriptReferences(uint sizeSkippedExon, uint sizeExons=200, uint transcriptNumber=2, uint exonNumber=8){
	vector<string> result;
	vector<string> exonList;
	// generate sequences for regular exons
	for(uint i(0); i < exonNumber; ++i){
		exonList.push_back(randomSequence(sizeExons));
	}
	// skipped exon
	exonList.push_back(randomSequence(sizeSkippedExon));
	// transcripts
	for(uint i(0); i < transcriptNumber; ++i){
		string transcript;
		transcript.reserve((exonNumber+1)*sizeExons);
		if (i == 0){  // inclusion isoform
			uint ii(0);
			while (ii < exonList.size()){
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






void generateReads(uint inclusionAbundance, vector<string>& transcripts, string& outSuffix, uint srCoverage=10, uint lrCoverage=100, uint errorRate=13, uint srLength=150){
	ofstream outLR("simulatedLR" + outSuffix + ".fa");
	ofstream outSR("simulatedSR" + outSuffix + ".fa");
	ofstream outRef("refSequences" + outSuffix + ".fa");
	ofstream outPerfectIncl("perfect_reads_inclusion" + outSuffix + ".fa");
	ofstream outPerfectExcl("perfect_reads_exclusion" + outSuffix + ".fa");
	ofstream outUncoIncl("uncorrected_reads_inclusion" + outSuffix + ".fa");
	ofstream outUncoExcl("uncorrected_reads_exclusion" + outSuffix + ".fa");
	string refSequence;
	uint numberSRInclusion(transcripts[0].size() / srLength * inclusionAbundance / 100 * srCoverage);
	uint numberSRExclusion(transcripts[1].size() / srLength * (100-inclusionAbundance) / 100 * srCoverage);
	uint numberLRInclusion(inclusionAbundance * lrCoverage / 100 );
	uint numberLRExclusion((100-inclusionAbundance) * lrCoverage / 100 );
	//~ uint numberSRInclusion(5);
	//~ uint numberSRExclusion(6);
	//~ uint numberLRInclusion(900);
	//~ uint numberLRExclusion((1000);
	cout << numberLRInclusion<< " " << numberLRExclusion << endl;
	cout << numberSRInclusion<< " " << numberSRExclusion << endl;
	string read;
	// LONG READS
	for (uint i(0); i < numberLRInclusion; ++i){
		read = mutateSequence(transcripts[0], errorRate);
		outLR << ">inclusion_" << i <<  endl << read << endl;
		outPerfectIncl << ">" << i << endl <<  transcripts[0] << endl;
		outUncoIncl << ">" << i << endl << read << endl;
	}
	for (uint i(0); i < numberLRExclusion; ++i){
		read = mutateSequence(transcripts[1], errorRate);
		outLR << ">exclusion_" << i << endl <<  read << endl;
		outPerfectExcl << ">" << i << endl << transcripts[1] << endl;
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
	outRef << ">inclusion" << endl << transcripts[0] << endl;
	outRef << ">exclusion" << endl << transcripts[1] << endl;
}






int main(int argc, char ** argv){
	if (argc > 2){
		uint sizeSkippedExon(stoi(argv[1]));
		uint inclusionAbundance(stoi(argv[2]));
		string outSuffix("");
		uint srCoverage(100), lrCoverage(10), errorRate(13);
		if (argc > 3){
			outSuffix = argv[3];
			if (argc > 5){
				srCoverage = stoi(argv[4]);
				lrCoverage = stoi(argv[5]);
				if (argc > 6){
					errorRate = stoi(argv[6]);
				}
			}
			
		}
		srand (time(NULL));
		auto startChrono = chrono::system_clock::now();
		// generate reference sequences for alternative transcripts
		vector<string> transcripts(generateAlternativeTranscriptReferences(sizeSkippedExon));
		// generate reads (long and short)
		generateReads(inclusionAbundance, transcripts, outSuffix, srCoverage, lrCoverage, errorRate);
		auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
		cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	} else {
		cout << "usage: ./ES_simulation <size skipped exon> <relative abudance inclusion isoform (%)> <error rate> " << endl;
	}
	return 0;
}
