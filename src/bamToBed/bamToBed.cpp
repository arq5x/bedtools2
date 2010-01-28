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
#define PROGRAM_NAME "bamToBed"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);
void PrintBed(const BamAlignment &bam, const RefVector &refs, bool useEditDistance);
void PrintBedPE(const BamAlignment &bam,  const RefVector &refs, bool useEditDistance);
bool IsCorrectMappingForBEDPE (const BamAlignment &bam);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string bamFile;

	bool haveBam = false;
	bool writeBedPE = false;
	bool useEditDistance = false;
	
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
		else if(PARAMETER_CHECK("-bedpe", 6, parameterLength)) {
				writeBedPE = true;
		}
		else if(PARAMETER_CHECK("-ed", 3, parameterLength)) {
				useEditDistance = true;
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
	// make sure we have an input files
	if (writeBedPE && useEditDistance) {
		cerr << endl << "*****" << endl << "*****ERROR: Cannot use edit distance with BEDPE. " << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (!showHelp) {

		// open the BAM file
		BamReader reader;
		reader.Open(bamFile);

		// get header & reference information
		string header = reader.GetHeaderText();
		RefVector refs = reader.GetReferenceData();
		
		BamAlignment bam;
		while (reader.GetNextAlignment(bam)) {	
			if (!writeBedPE) {
				if (bam.IsMapped()) {
					PrintBed(bam, refs, useEditDistance);
				}
			}
			else {
				if ( ! bam.IsPaired() ) {	// ensure that the BAM file is paired.
					cerr << "Encountered an unpaired alignment.  Are you sure this is a paired BAM file?  Exiting." << endl;
					exit(1);
				}
				else if (IsCorrectMappingForBEDPE(bam)) {	// only report one of the two mappings for a pair
					PrintBedPE(bam, refs, useEditDistance);
				}
			}	
		}
		reader.Close();	
	}	
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {

	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Converts BAM alignments to BED6 or BEDPE format." << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;

	cerr << "Options: " << endl;
	
	cerr << "\t-bedpe\t"	<< "Write BEDPE format." << endl << endl;
	
	cerr << "\t-ed\t"		<< "Use BAM edit distance (NM tag) for score." << endl;
	cerr 					<< "\t\tDefault is to use mapping quality." << endl;
	cerr 					<< "\t\tNot available for BEDPE format." << endl << endl;
	

	// end the program here
	exit(1);

}


void PrintBed(const BamAlignment &bam,  const RefVector &refs, bool useEditDistance) {

	string strand = "+"; 
	if (bam.IsReverseStrand()) strand = "-";
	string name = bam.Name;
	if (bam.IsFirstMate()) name += "/1";
	if (bam.IsSecondMate()) name += "/2";
	
	if (useEditDistance == false) {
		printf("%s\t%d\t%d\t\%s\t%d\t%s\n", refs.at(bam.RefID).RefName.c_str(), bam.Position,
									  bam.Position + bam.Length, name.c_str(),
									  bam.MapQuality, strand.c_str());
	}
	else {
		uint8_t editDistance;
		if (bam.GetEditDistance(editDistance)) {
			printf("%s\t%d\t%d\t\%s\t%u\t%s\n", refs.at(bam.RefID).RefName.c_str(), bam.Position,
										  bam.Position + bam.Length, name.c_str(),
										  editDistance, strand.c_str());
		}
		else {
			cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
			exit(1);
		}
	}
}


void PrintBedPE(const BamAlignment &bam,  const RefVector &refs, bool useEditDistance) {

	string chrom1, chrom2, strand1, strand2;
	int start1, start2, end1, end2;
	start1 = start2 = end1 = end2 = -1;
	chrom1 = chrom2 = strand1 = strand2 = ".";
	
	// end 1
	if (bam.IsMapped()) {
		chrom1 = refs.at(bam.RefID).RefName;
		start1 = bam.Position;
		end1 = bam.Position + bam.Length;
		strand1 = "+";
		if (bam.IsReverseStrand()) strand1 = "-";	
	}

	// end 2
	if (bam.IsMateMapped()) {
		chrom2 = refs.at(bam.MateRefID).RefName;
		start2 = bam.MatePosition;
		end2 = bam.MatePosition + bam.Length;
		strand2 = "+";
		if (bam.IsMateReverseStrand()) strand2 = "-";	
	}		

	printf("%s\t%d\t%d\t\%s\t%d\t%d\t%s\t%d\t%s\t%s\n", 
			chrom1.c_str(), start1, end1, chrom2.c_str(), start2, end2, 
			bam.Name.c_str(), bam.MapQuality, strand1.c_str(), strand2.c_str());

}



bool IsCorrectMappingForBEDPE (const BamAlignment &bam) {

	if ( (bam.RefID == bam.MateRefID) && (bam.InsertSize > 0) ) {
		return true;
	}
	else if ( (bam.RefID == bam.MateRefID) && (bam.InsertSize == 0) && bam.IsFirstMate() ) {
		return true;
	}
	else if ( (bam.RefID != bam.MateRefID) && bam.IsFirstMate() ) {
		return true;
	}
	else return false;
}
