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

	_bedFile = bedFile;
	_genomeFile = genomeFile;
	_forceStrand = forceStrand;
	
	_leftSlop = leftSlop;
	_rightSlop = rightSlop;
	
	_bed    = new BedFile(bedFile);
	_genome = new GenomeFile(genomeFile);
	
	SlopBed();
}


BedSlop::~BedSlop(void) {

}


void BedSlop::SlopBed() {
	
	int lineNum = 0;
	BED bedEntry, nullBed;     // used to store the current BED line from the BED file.

	_bed->Open();
	while (_bed->GetNextBed(bedEntry, lineNum)) {
		AddSlop(bedEntry);
		_bed->reportBedNewLine(bedEntry);
		bedEntry = nullBed;					
	}
	_bed->Close();
}


void BedSlop::AddSlop(BED &bed) {

	// special handling if the BED entry is on the negative
	// strand and the user cares about strandedness.
	int chromSize = _genome->getChromSize(bed.chrom);
	
	if ( (_forceStrand) && (bed.strand == "-") ) {
		// inspect the start
		if ((bed.start - _rightSlop) > 0) bed.start -= _rightSlop;
		else bed.start = 0;

		// inspect the start		
		if ((bed.end + _leftSlop) <= chromSize) bed.end += _leftSlop;
		else bed.end = chromSize;
	}
	else {		
		// inspect the start
		if ((bed.start - _leftSlop) > 0) bed.start -= _leftSlop;
		else bed.start = 0;
		
		// inspect the end
		if ((bed.end + _rightSlop) <= chromSize) bed.end += _rightSlop;
		else bed.end = chromSize;
	}
}


