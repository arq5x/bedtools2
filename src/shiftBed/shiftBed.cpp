/*****************************************************************************
  shiftBed.cpp

  (c) 2016 - David Richardson
  European Molecular Biology Laboratory, European Bioinformatics Institute
  davidr@ebi.ac.uk

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "shiftBed.h"

BedShift::BedShift(string &bedFile, string &genomeFile, 
            float shiftMinus, float shiftPlus, bool fractional, bool printHeader) {

    _bedFile     = bedFile;
    _genomeFile  = genomeFile;
    _shiftMinus  = shiftMinus;
    _shiftPlus   = shiftPlus;
    _fractional  = fractional;
    _printHeader = printHeader;

    _bed    = new BedFile(bedFile);
    _genome = new GenomeFile(genomeFile);

    // get going, shift it around.
    ShiftBed();
}


BedShift::~BedShift(void) {

}


void BedShift::ShiftBed() {

    BED bedEntry;     // used to store the current BED line from the BED file.
    float m, p;

    _bed->Open();
    // report header first if asked.
    if (_printHeader == true) {
        _bed->PrintHeader();
    }        
    while (_bed->GetNextBed(bedEntry)) {    
        if (_bed->_status == BED_VALID) {
            if (_fractional == false) {
                AddShift(bedEntry);
            }
            else {
	        m = _shiftMinus;	
                _shiftMinus  = _shiftMinus * (float)bedEntry.size();
	        p = _shiftPlus;	
                _shiftPlus = _shiftPlus * (float)bedEntry.size();
                AddShift(bedEntry);
	        _shiftMinus = m;
	        _shiftPlus = p;
            }
            _bed->reportBedNewLine(bedEntry);
        }
    }
    _bed->Close();
}


void BedShift::AddShift(BED &bed) {

    CHRPOS chromSize = (CHRPOS)_genome->getChromSize(bed.chrom);

	int start,end;

	if (bed.strand == "-"){
		start = bed.start + (int)_shiftMinus;
		end = bed.end + (int)_shiftMinus;
	}
	else {
		start = bed.start + (int)_shiftPlus;
		end = bed.end + (int)_shiftPlus;
	}
	
	// has the entry run off the end of the genome?
	if (start < 0) {
		start = 0;
	}
	if (end < 1) {
		end = 1;
	}
	if (start > (chromSize-1)){
		start = (chromSize-1);
	}
	if (end > chromSize){
		end = chromSize;
	}
	
	bed.start = start;
	bed.end = end;

}


