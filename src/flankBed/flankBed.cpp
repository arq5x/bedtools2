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


BedFlank::BedFlank(string &bedFile, string &genomeFile, bool forceStrand, 
                   float leftFlank, float rightFlank, bool fractional,
                   bool printHeader) 
{

    _bedFile      = bedFile;
    _genomeFile   = genomeFile;
    _forceStrand  = forceStrand;
    _leftFlank    = leftFlank;
    _rightFlank   = rightFlank;
    _fractional   = fractional;
    _printHeader  = printHeader;

    _bed    = new BedFile(bedFile);
    _genome = new GenomeFile(genomeFile);

    // get going, slop it up.
    FlankBed();
}


BedFlank::~BedFlank(void) {

}


void BedFlank::FlankBed() {

    BED bedEntry;     // used to store the current BED line from the BED file.

    _bed->Open();
    // report A's header first if asked.
    if (_printHeader == true) {
        _bed->PrintHeader();
    }
        
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            int leftFlank  = _leftFlank;
            int rightFlank = _rightFlank;            
            if (_fractional == true) {
                leftFlank  = (int) (_leftFlank  * bedEntry.size());
                rightFlank = (int) (_rightFlank * bedEntry.size());
            }
            
            if ((_forceStrand == false) || (bedEntry.strand == "+"))
            {
                AddFlank(bedEntry,  leftFlank, rightFlank);
            }
            else if ((_forceStrand == true) && (bedEntry.strand == "-" ))
            {
                AddStrandedFlank(bedEntry,  leftFlank, rightFlank);                    
            }
        }
    }
    _bed->Close();
}


void BedFlank::AddFlank(BED &bed, int leftFlank, int rightFlank) {

    int chromSize = _genome->getChromSize(bed.chrom);
    if (chromSize == -1) {
        cerr << "ERROR: chrom \"" << bed.chrom << "\" not found in genome file. Exiting." << endl;
        exit(1);
    }

    // init. our left and right flanks to the original BED entry.
    // we'll create the flanks from these coordinates.
    BED left  = bed;
    BED right = bed;
    // split of feature bed to be similar to BedFlank::AddStrandedFlank
    BED feature = bed;
    
    left.zeroLength = false;
    right.zeroLength = false;
    
    // make the left flank (if necessary)
    // do not report a left flank if the left feature start is already
    // the very first base of the scaffold
    if ((leftFlank > 0) && (feature.start != 0)) {
        if ( (static_cast<int>(feature.start) - leftFlank) > 0)
        {
            left.end    = feature.start;
            left.start -= leftFlank;
        }
        else 
        {
            left.end    = feature.start;
            left.start  = 0;
        }
        // report the left flank
        _bed->reportBedNewLine(left);
    }
    
    // make the right flank (if necessary)
    // Do not flank if the right flank of the feature (which is feature.end)
    // is already last base of scaffold
    if (rightFlank > 0 && (feature.end < chromSize) ) {
        if ( (static_cast<int>(feature.end) + static_cast<int>(rightFlank+1))
             <= static_cast<int>(chromSize)) 
        {
            right.start    = feature.end;
            right.end     += (rightFlank);
        }
        else {
            right.start    = feature.end;
            right.end      = chromSize;
        }
        // report the right flank
        _bed->reportBedNewLine(right);
    }    
}


void BedFlank::AddStrandedFlank(BED &bed, int leftFlank, int rightFlank) {

    int chromSize = _genome->getChromSize(bed.chrom);
    if (chromSize == -1) {
        cerr << "ERROR: chrom \"" << bed.chrom << "\" not found in genome file. Exiting." << endl;
        exit(1);
    }

    // init. our left and right flanks to the original BED entry.
    BED left  = bed;
    BED right = bed;
    // separate feature bed file because when calling both left and right flanks otherwise gives conflicts
    BED feature = bed;
    
    // make the right flank (if necessary)
    // Do not flank if the right flank of the feature (which is feature.start)
    // since in reverse strand is first base of scaffold
    if (rightFlank > 0 && (feature.start != 0)) {
        if ( (static_cast<int>(feature.start) - rightFlank) > 0)
        {
            left.end    = feature.start;
            left.start -= rightFlank;
        }
        else 
        {
            left.end    = feature.start;
            left.start  = 0;
        }
        // report the left flank
        _bed->reportBedNewLine(left);
    }
    
    // make the left flank (if necessary)
    // left.start is on the right side of the feature because feature is on reverse strand
    // do not report a left flank if the feature.end is already last base of the scaffold
    // coordinates here are BED so end is exclusive
    // in BED coords, this will create a zero length feature with flank.start == flank.end
    // however when converting however to GFF, the start will get +1
    // and so the feature start will be 1 higher than the end
    // resulting in errors in further bed tools
    // you would create a left flank that has as left.start the feature.end coordinate
    // so end needs to be smaller than chromsize to have at least 1 flanking base
    if (leftFlank > 0 && (feature.end < chromSize) ) {
        if ( (static_cast<int>(feature.end) + leftFlank) <= static_cast<int>(chromSize))
        {
            right.start    = feature.end;
            right.end     += leftFlank;
        }
        else {
            right.start    = feature.end;
            right.end      = chromSize;
        }
        // report the right flank
        _bed->reportBedNewLine(right);
    }   
}


