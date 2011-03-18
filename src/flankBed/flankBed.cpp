/*****************************************************************************
  flankBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "flankBed.h"


BedFlank::BedFlank(string &bedFile, string &genomeFile, bool forceStrand, float leftFlank, float rightFlank, bool fractional) {

    _bedFile      = bedFile;
    _genomeFile   = genomeFile;
    _forceStrand  = forceStrand;
    _leftFlank    = leftFlank;
    _rightFlank   = rightFlank;
    _fractional   = fractional; 

    _bed    = new BedFile(bedFile);
    _genome = new GenomeFile(genomeFile);

    // get going, slop it up.
    FlankBed();
}


BedFlank::~BedFlank(void) {

}


void BedFlank::FlankBed() {

    int lineNum = 0;
    BED bedEntry, nullBed;     // used to store the current BED line from the BED file.
    BedLineStatus bedStatus;

    _bed->Open();
    bedStatus = _bed->GetNextBed(bedEntry, lineNum);
    while (bedStatus != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            if (_fractional == false) {
                if (_forceStrand == false) {
                    AddFlank(bedEntry,  (int) _leftFlank, (int) _rightFlank);
                }
                else {
                    AddStrandedFlank(bedEntry,  (int) _leftFlank, (int) _rightFlank);                    
                }
            }
            else {
                int leftFlank  = (int) (_leftFlank  * bedEntry.size());
                int rightFlank = (int) (_rightFlank * bedEntry.size());
                if (_forceStrand == false) {
                    AddFlank(bedEntry, leftFlank, rightFlank);
                }
                else {
                    AddStrandedFlank(bedEntry, leftFlank, rightFlank);                    
                }
            }
            bedEntry = nullBed;
        }
        bedStatus = _bed->GetNextBed(bedEntry, lineNum);
    }
    _bed->Close();
}


void BedFlank::AddFlank(BED &bed, int leftFlank, int rightFlank) {

    CHRPOS chromSize = _genome->getChromSize(bed.chrom);

    // init. our left and right flanks to the original BED entry.
    // we'll create the flanks from these coordinates.
    BED left  = bed;
    BED right = bed;
    
    // make the left flank (if necessary)
    if (leftFlank > 0) {
        if ( (static_cast<int>(left.start) - leftFlank) > 0) 
        {
            left.end    = left.start;
            left.start -= leftFlank;
        }
        else 
        {
            left.end    = left.start;
            left.start  = 0;
        }
        // report the left flank
        _bed->reportBedNewLine(left);
    }
    
    // make the left flank (if necessary)
    if (rightFlank > 0) {
        if ( (static_cast<int>(right.end) + (rightFlank+1)) <= static_cast<int>(chromSize)) 
        {
            right.start    = right.end;
            right.end     += (rightFlank);
        }
        else {
            right.start    = right.end;
            right.end     += chromSize;
        }
        // report the right flank
        _bed->reportBedNewLine(right);
    }    
}


void BedFlank::AddStrandedFlank(BED &bed, int leftFlank, int rightFlank) {

    CHRPOS chromSize = _genome->getChromSize(bed.chrom);

    // init. our left and right flanks to the original BED entry.
    // we'll create the flanks from these coordinates.
    BED left  = bed;
    BED right = bed;
    
    // make the left flank (if necessary)
    if (rightFlank > 0) {
        if ( (static_cast<int>(left.start) - rightFlank) > 0) 
        {
            left.end    = left.start;
            left.start -= rightFlank;
        }
        else 
        {
            left.end    = left.start;
            left.start  = 0;
        }
        // report the left flank
        _bed->reportBedNewLine(left);
    }
    
    // make the left flank (if necessary)
    if (leftFlank > 0) {
        if ( (static_cast<int>(right.end) + (leftFlank+1)) <= static_cast<int>(chromSize)) 
        {
            right.start    = right.end;
            right.end     += (leftFlank);
        }
        else {
            right.start    = right.end;
            right.end      = chromSize;
        }
        // report the right flank
        _bed->reportBedNewLine(right);
    }   
}


