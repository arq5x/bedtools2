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

/*
	Constructor
*/
BedSlop::BedSlop(string &bedFile, string &genomeFile, bool &forceStrand, int &leftSlop, int &rightSlop) {

	this->bedFile = bedFile;
	this->genomeFile = genomeFile;
	this->forceStrand = forceStrand;
	
	this->leftSlop = leftSlop;
	this->rightSlop = rightSlop;
	
	this->bed = new BedFile(bedFile);	
}



/*
	Destructor
*/
BedSlop::~BedSlop(void) {

}


void BedSlop::SlopBed(istream &bedInput) {
	
	BED bedEntry;     // used to store the current BED line from the BED file.
	int lineNum = 0;
	string bedLine;	  // used to store the current (unparsed) line from the BED file.
		
	while (getline(bedInput, bedLine)) {
		
		vector<string> bedFields;
		Tokenize(bedLine,bedFields);
		lineNum++;
		
		/*
		   if a valid BED entry, add requested "slop" and print out
		   the adjusted BED entry.
		*/
		if (this->bed->parseLine(bedEntry, bedFields, lineNum)) {
			AddSlop(bedEntry);
			bed->reportBedNewLine(bedEntry);			
		}
	}
}


void BedSlop::AddSlop(BED &bed) {

	/*
	   special handling if the BED entry is on the negative
	   strand and the user cares about strandedness.
	*/
	if ((this->forceStrand) &&  (bed.strand == "-")) {
		// inspect the start
		if ((bed.start - rightSlop) > 0) bed.start -= rightSlop;
		else bed.start = 0;

		// inspect the start		
		if ((bed.end + leftSlop) <= this->chromSizes[bed.chrom]) bed.end += leftSlop;
		else bed.end = this->chromSizes[bed.chrom];
	}
	else {		
		// inspect the start
		if ((bed.start - leftSlop) > 0) bed.start -= leftSlop;
		else bed.start = 0;
		
		// inspect the end
		if ((bed.end + rightSlop) <= this->chromSizes[bed.chrom]) bed.end += rightSlop;
		else bed.end = this->chromSizes[bed.chrom];
	}
}


void BedSlop::DetermineBedInput() {


	/* open the GENOME file for reading.
	   if successful, load each chrom and it's size into
	   the "chromSize" map.  also compute the total size of the genome
	   and store in "genomeSize"
	*/
	ifstream genome(this->genomeFile.c_str(), ios::in);
	if ( !genome ) {
		cerr << "Error: The requested genome file (" <<this->genomeFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	else {
		string chrom;
		unsigned int size;
		while (genome >> chrom >> size) {
			if (chrom.size() > 0 && size > 0) {
				this->chromSizes[chrom] = size;
			}
		}
	}

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
