/*****************************************************************************
  bamToBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "version.h"
#include "sequenceUtils.h"
#include "BamReader.h"
#include "BamWriter.h"
#include "BamAux.h"
using namespace BamTools;

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bamFillMateSeq"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);
void GetSeqsAndQuals(const vector<BamAlignment> &alignments,
                     string &seq1,
                     string &qual1,
                     string &seq,
					 string &qual2);
void FillMateSeqForReadBlock (const vector<BamAlignment> &alignments, 
                             const RefVector &refs,
                             BamWriter &writer,
		                     const string &seq1,
		                     const string &qual1,
		                     const string &seq,
							 const string &qual2);
int strnum_cmp(const char *a, const char *b);

	
int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;
	bool haveBam  = false;
	
	// input files
	string bamFile;
			
	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
		(PARAMETER_CHECK("--help", 5, parameterLength))) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();

	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-i", 2, parameterLength)) {
			if ((i+1) < argc) {
				haveBam = true;
				bamFile = argv[i + 1];
				i++;
			}
		}
		else {
			cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have an input files
	if (!haveBam ) {
		cerr << endl << "*****" << endl << "*****ERROR: Need -i (BAM) file. " << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (!showHelp) {

		// open the BAM file
		BamReader reader;
		reader.Open(bamFile);
		BamWriter writer;

		// get header & reference information
		string header = reader.GetHeaderText();
		RefVector refs = reader.GetReferenceData();

		// open our BAM writer
		writer.Open("stdout", header, refs);
		
		// track the previous and current sequence
		// names so that we can identify blocks of
		// alignments for a given read ID.
		string prevName = "";
		string currName = "";
		
		vector<BamAlignment> alignments;		// vector of BAM alignments for a given ID in a BAM file.
		alignments.reserve(100);
		
		BamAlignment bam;
		string seq1, qual1, seq2, qual2;
		// get each set of alignments for each pair.
		while (reader.GetNextAlignment(bam)) {
		
			currName = bam.Name;
			//printf("%s\t%s\n", currName.c_str(), prevName.c_str());		
			if ( (strcmp(currName.c_str(), prevName.c_str()) != 0) && (prevName != "") ) {
				//printf("yep\n");
				// block change.  enforce expected sort order.
				if (strnum_cmp(currName.c_str(), prevName.c_str()) == -1) {
					cerr << "Error: The BAM file (" << bamFile << ") is not sorted by sequence name. Exiting!" << endl;
					exit (1);
				}
				else {
					// process the current block of alignments
					GetSeqsAndQuals(alignments, seq1, qual1, seq2, qual2);
					FillMateSeqForReadBlock(alignments, refs, writer, seq1, qual1, seq2, qual2);
										
					// reset the alignment vector for the next block
					// and add the first alignment in the next block
					alignments.clear();
					alignments.push_back(bam);
				}	
			}
			else {
				alignments.push_back(bam);
			}
			prevName = currName;
		}
		// process the last block
		GetSeqsAndQuals(alignments, seq1, qual1, seq2, qual2);
		FillMateSeqForReadBlock(alignments, refs, writer, seq1, qual1, seq2, qual2);

		// close up
		reader.Close();	
		writer.Close();
	}	
	else {
		ShowHelp();
	}
}


void GetSeqsAndQuals(const vector<BamAlignment> &alignments,
                     string &seq1,
                     string &qual1,
                     string &seq2,
                     string &qual2) {

	bool haveSeq1, haveSeq2, haveQual1, haveQual2;
	haveSeq1 = haveSeq2 = haveQual1 = haveQual2 = false;

	// extract the sequences for the first and second mates

	vector<BamAlignment>::const_iterator alItr = alignments.begin();
	vector<BamAlignment>::const_iterator alEnd = alignments.end();
	while ( (alItr != alEnd) && (haveSeq1 == false || haveSeq2 == false || 
	                             haveQual1 == false || haveQual2 == false)) {

		if ( alItr->IsFirstMate() ) {
			seq1  = alItr->QueryBases;
			qual1 = alItr->Qualities;
			if ( alItr->IsReverseStrand() ) {
				reverseComplement(seq1);
				reverseSequence(qual1);
			}
			haveSeq1  = true;
			haveQual1 = true;
		}
		else if ( alItr->IsSecondMate() ) {
			seq2  = alItr->QueryBases;
			qual2 = alItr->Qualities;
			if ( alItr->IsReverseStrand() ) {
				reverseComplement(seq2);
				reverseSequence(qual2);
			}
			haveSeq2  = true;
			haveQual2 = true; 
		}
		++alItr;
	}

	//cerr << alignments.size() << endl;
	cout << alignments[0].Qualities << endl;
	//seq1  = alignments[0].QueryBases;
	//qual1 = alignments[0].Qualities;
	//seq2  = alignments[1].QueryBases;
	//qual2 = alignments[1].Qualities;
	seq1 = "seq1\0";
	seq2 = "seq2\0";
	qual1 = "88###########################################################################88";
	qual2 = "77###########################################################################77";
	
	// IT MUST have to do with modifying the alignment itself.
}


void FillMateSeqForReadBlock (const vector<BamAlignment> &alignments, 
                             const RefVector &refs,
                             BamWriter &writer,
		                     const string &seq1,
		                     const string &qual1,
			                 const string &seq2,
		                     const string &qual2) {
	// now update each alignment with the appropriate 
	// sequence for it's mate
	cerr << "call\n";
	string joe = qual1;
	vector<BamAlignment>::const_iterator aItr = alignments.begin();
	vector<BamAlignment>::const_iterator aEnd = alignments.end();
	int i = 0;
	for (; aItr != aEnd; ++aItr) {
		
		cerr << i << " " << aItr->Name << " [" << qual1 << "]" << "\t[" << qual2 << "]" << "\t" << joe << endl;
		//if ( aItr->IsFirstMate() ) {
			//aItr->AddBamTag("R2", "Z", seq2);
			//aItr->AddBamTag("Q2", "Z", qual2);			
		//}
		//else if ( aItr->IsSecondMate() ) {
			//aItr->AddBamTag("R2", "Z", seq1);
			//aItr->AddBamTag("Q2", "Z", qual1);			
		//}
		//BamAlignment poop = *aItr;
		writer.SaveAlignment(*aItr);
		cerr << i << " " <<  aItr->Name << " [" << qual1 << "]" << "\t[" << qual2 << "]" <<  "\t" << joe << endl << endl;
		i++;
	}
}


// samtools sort comparison function from bam_sort.c
// used for sorting by query name.
int strnum_cmp(const char *a, const char *b) {
	char *pa, *pb;
	pa = (char*)a; pb = (char*)b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
		} else {
			if (*pa != *pb) break;
			++pa; ++pb;
		}
	}
	if (*pa == *pb)
		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}


void ShowHelp(void) {

	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Updates the R2 and Q2 BAM tag to include the mate's seq/qual." << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;

//	cerr << "Options: " << endl;
	
//	cerr << "\t-bedpe\t"	<< "Write BEDPE format." << endl << endl;

	cerr << "Notes: " << endl;

	cerr << "\tInput BAM must be sorted by read name" << endl << endl;

	// end the program here
	exit(1);
}

