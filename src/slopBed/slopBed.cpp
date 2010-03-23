/*****************************************************************************
  slopBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "slopBed.h"


BedSlop::BedSlop(string &bedFile, string &genomeFile, bool &forceStrand, int &leftSlop, int &rightSlop) {

	this->bedFile = bedFile;
	this->genomeFile = genomeFile;
	this->forceStrand = forceStrand;
	
	this->leftSlop = leftSlop;
	this->rightSlop = rightSlop;
	
	this->bed    = new BedFile(bedFile);
	this->genome = new GenomeFile(genomeFile);	
}


BedSlop::~BedSlop(void) {

}


void BedSlop::SlopBed(istream &bedInput) {
	
	int lineNum = 0;
	string bedLine;	  // used to store the current (unparsed) line from the BED file.
		
	while (getline(bedInput, bedLine)) {
		
		vector<string> bedFields;
		Tokenize(bedLine,bedFields);
		lineNum++;
		
		BED bedEntry;     // used to store the current BED line from the BED file.
		if (this->bed->parseLine(bedEntry, bedFields, lineNum)) {
			AddSlop(bedEntry);
			bed->reportBedNewLine(bedEntry);			
		}
	}
}


void BedSlop::AddSlop(BED &bed) {

	// special handling if the BED entry is on the negative
	// strand and the user cares about strandedness.
	int chromSize = genome->getChromSize(bed.chrom);
	
	if ( (this->forceStrand) && (bed.strand == "-") ) {
		// inspect the start
		if ((bed.start - rightSlop) > 0) bed.start -= rightSlop;
		else bed.start = 0;

		// inspect the start		
		if ((bed.end + leftSlop) <= chromSize) bed.end += leftSlop;
		else bed.end = chromSize;
	}
	else {		
		// inspect the start
		if ((bed.start - leftSlop) > 0) bed.start -= leftSlop;
		else bed.start = 0;
		
		// inspect the end
		if ((bed.end + rightSlop) <= chromSize) bed.end += rightSlop;
		else bed.end = chromSize;
	}
}


void BedSlop::DetermineBedInput() {

	if (this->bedFile != "stdin") {   // process a file
		ifstream beds(this->bedFile.c_str(), ios::in);
		if ( !beds ) {
			cerr << "Error: The requested bed file (" << this->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		SlopBed(beds);
	}
	else {   // process stdin
		SlopBed(cin);		
	}
}
