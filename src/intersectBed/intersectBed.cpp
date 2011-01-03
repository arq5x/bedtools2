/*****************************************************************************
  intersectBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "intersectBed.h"

/************************************
Helper functions
************************************/
bool BedIntersect::processHits(const BED &a, const vector<BED> &hits, bool printable) {

    // how many overlaps are there b/w the bed and the set of hits?
    int s, e, overlapBases;
    int  numOverlaps = 0;
    bool hitsFound   = false;
    int aLength      = (a.end - a.start);   // the length of a in b.p.

    // loop through the hits and report those that meet the user's criteria
    vector<BED>::const_iterator h       = hits.begin();
    vector<BED>::const_iterator hitsEnd = hits.end();
    for (; h != hitsEnd; ++h) {
        s            = max(a.start, h->start);
        e            = min(a.end, h->end);
        overlapBases = (e - s);             // the number of overlapping bases b/w a and b

        // is there enough overlap relative to the user's request? (default ~ 1bp)
        if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) {
            // Report the hit if the user doesn't care about reciprocal overlap between A and B.
            if (_reciprocal == false) {
                hitsFound = true;
                numOverlaps++;
                if (printable == true)
                    ReportOverlapDetail(overlapBases, a, *h, s, e);
            }
            // we require there to be sufficient __reciprocal__ overlap
            else {
                int bLength    = (h->end - h->start);
                float bOverlap = ( (float) overlapBases / (float) bLength );
                if (bOverlap >= _overlapFraction) {
                    hitsFound = true;
                    numOverlaps++;
                    if (printable == true)
                        ReportOverlapDetail(overlapBases, a, *h, s, e);
                }
            }
        }
    }
    // report the summary of the overlaps if requested.
    ReportOverlapSummary(a, numOverlaps);
    // were hits found for this BED feature?
    return hitsFound;
}


/*
    Constructor
*/
BedIntersect::BedIntersect(string bedAFile, string bedBFile, bool anyHit,
                           bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                           float overlapFraction, bool noHit, bool writeCount, bool forceStrand,
                           bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput, bool isUncompressedBam) {

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _anyHit              = anyHit;
    _noHit               = noHit;
    _writeA              = writeA;
    _writeB              = writeB;
    _writeOverlap        = writeOverlap;
    _writeAllOverlap     = writeAllOverlap;
    _writeCount          = writeCount;
    _overlapFraction     = overlapFraction;
    _forceStrand         = forceStrand;
    _reciprocal          = reciprocal;
    _obeySplits          = obeySplits;
    _bamInput            = bamInput;
    _bamOutput           = bamOutput;
    _isUncompressedBam   = isUncompressedBam;

    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);

    if (_bamInput == false)
        IntersectBed();
    else
        IntersectBam(bedAFile);
}


/*
    Destructor
*/
BedIntersect::~BedIntersect(void) {
}


bool BedIntersect::FindOverlaps(const BED &a, vector<BED> &hits) {

    bool hitsFound = false;

    // should we print each overlap, or does the user want summary information?
    bool printable = true;
    if (_anyHit || _noHit || _writeCount)
        printable = false;

    // collect and report the sufficient hits
    _bedB->FindOverlapsPerBin(a.chrom, a.start, a.end, a.strand, hits, _forceStrand);
    hitsFound = processHits(a, hits, printable);

    return hitsFound;
}


void BedIntersect::ReportOverlapDetail(const int &overlapBases, const BED &a, const BED &b,
                                       const CHRPOS &s, const CHRPOS &e) {
    // default. simple intersection only
    if (_writeA == false && _writeB == false && _writeOverlap == false) {
        _bedA->reportBedRangeNewLine(a,s,e);
    }
    //  -wa -wbwrite the original A and B
    else if (_writeA == true && _writeB == true) {
        _bedA->reportBedTab(a);
        _bedB->reportBedNewLine(b);
    }
    // -wa write just the original A
    else if (_writeA == true) {
        _bedA->reportBedNewLine(a);
    }
    // -wb write the intersected portion of A and the original B
    else if (_writeB == true) {
        _bedA->reportBedRangeTab(a,s,e);
        _bedB->reportBedNewLine(b);
    }
    // -wo write the original A and B plus the no. of overlapping bases.
    else if (_writeOverlap == true) {
        _bedA->reportBedTab(a);
        _bedB->reportBedTab(b);
        printf("%d\n", overlapBases);
    }
}


void BedIntersect::ReportOverlapSummary(const BED &a, const int &numOverlapsFound) {
    // -u  just report the fact that there was >= 1 overlaps
    if (_anyHit && (numOverlapsFound >= 1)) {
        _bedA->reportBedNewLine(a);
    }
    // -c  report the total number of features overlapped in B
    else if (_writeCount) {
        _bedA->reportBedTab(a);
        printf("%d\n", numOverlapsFound);
    }
    // -v  report iff there were no overlaps
    else if (_noHit && (numOverlapsFound == 0)) {
        _bedA->reportBedNewLine(a);
    }
    // -wao the user wants to force the reporting of 0 overlap
    else if (_writeAllOverlap && (numOverlapsFound == 0)) {
        _bedA->reportBedTab(a);
        _bedB->reportNullBedTab();
        printf("0\n");
    }
}


bool BedIntersect::FindOneOrMoreOverlap(const BED &a) {
    bool overlapsFound;
    if (_reciprocal == false) {
        overlapsFound = _bedB->FindOneOrMoreOverlapsPerBin(a.chrom, a.start, a.end, a.strand,
                                                          _forceStrand, _overlapFraction);
    }
    else {
        overlapsFound = _bedB->FindOneOrMoreReciprocalOverlapsPerBin(a.chrom, a.start, a.end, a.strand,
                                                                    _forceStrand, _overlapFraction);
    }
    return overlapsFound;
}


void BedIntersect::IntersectBed() {

    // load the "B" file into a map in order to
    // compare each entry in A to it in search of overlaps.
    _bedB->loadBedFileIntoMap();

    int lineNum = 0;
    vector<BED> hits;
    hits.reserve(100);
    BED a, nullBed;
    BedLineStatus bedStatus;

    // open the "A" file, process each BED entry and searh for overlaps.
    _bedA->Open();
    while ((bedStatus = _bedA->GetNextBed(a, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            // treat the BED as a single "block"
            if (_obeySplits == false) {
                FindOverlaps(a, hits);
                hits.clear();
                a = nullBed;
            }
            // split the BED12 into blocks and look for overlaps in each discrete block
            else {
                bedVector bedBlocks;  // vec to store the discrete BED "blocks"
                splitBedIntoBlocks(a, lineNum, bedBlocks);

                vector<BED>::const_iterator bedItr  = bedBlocks.begin();
                vector<BED>::const_iterator bedEnd  = bedBlocks.end();
                for (; bedItr != bedEnd; ++bedItr) {
                    FindOverlaps(*bedItr, hits);
                    hits.clear();
                }
                a = nullBed;
            }
        }
    }
    _bedA->Close();
}


void BedIntersect::IntersectBam(string bamFile) {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedFileIntoMap();

    // open the BAM file
    BamReader reader;
    BamWriter writer;
    reader.Open(bamFile);

    // get header & reference information
    string header  = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // open a BAM output to stdout if we are writing BAM
    if (_bamOutput == true) {
        // open our BAM writer
        writer.Open("stdout", header, refs, _isUncompressedBam);
    }

    vector<BED> hits;
    // reserve some space
    hits.reserve(100);

    _bedA->bedType = 6;
    BamAlignment bam;
    // get each set of alignments for each pair.
    while (reader.GetNextAlignment(bam)) {

        if (bam.IsMapped()) {
            BED a;
            a.chrom = refs.at(bam.RefID).RefName;
            a.start = bam.Position;
            a.end   = bam.GetEndPosition(false);

            // build the name field from the BAM alignment.
            a.name = bam.Name;
            if (bam.IsFirstMate()) a.name += "/1";
            if (bam.IsSecondMate()) a.name += "/2";

            a.score  = ToString(bam.MapQuality);

            a.strand = "+";
            if (bam.IsReverseStrand()) a.strand = "-";

            if (_bamOutput == true) {
                bool overlapsFound = false;
                // treat the BAM alignment as a single "block"
                if (_obeySplits == false) {
                    overlapsFound = FindOneOrMoreOverlap(a);
                }
                // split the BAM alignment into discrete blocks and
                // look for overlaps only within each block.
                else {
                    bool overlapFoundForBlock;
                    bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
                    // we don't want to split on "D" ops, hence the "false"
                    getBamBlocks(bam, refs, bedBlocks, false);

                    vector<BED>::const_iterator bedItr  = bedBlocks.begin();
                    vector<BED>::const_iterator bedEnd  = bedBlocks.end();
                    for (; bedItr != bedEnd; ++bedItr) {
                        overlapFoundForBlock = FindOneOrMoreOverlap(a);
                        if (overlapFoundForBlock == true)
                            overlapsFound = true;
                    }
                }
                if (overlapsFound == true) {
                    if (_noHit == false)
                        writer.SaveAlignment(bam);
                }
                else {
                    if (_noHit == true) {
                        writer.SaveAlignment(bam);
                    }
                }
            }
            else {
                // treat the BAM alignment as a single BED "block"
                if (_obeySplits == false) {
                    FindOverlaps(a, hits);
                    hits.clear();
                }
                // split the BAM alignment into discrete BED blocks and
                // look for overlaps only within each block.
                else {
                    bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
                    getBamBlocks(bam, refs, bedBlocks, false);

                    vector<BED>::const_iterator bedItr  = bedBlocks.begin();
                    vector<BED>::const_iterator bedEnd  = bedBlocks.end();
                    for (; bedItr != bedEnd; ++bedItr) {
                        FindOverlaps(*bedItr, hits);
                        hits.clear();
                    }
                }
            }
        }
    }

    // close the relevant BAM files.
    reader.Close();
    if (_bamOutput == true) {
        writer.Close();
    }
}

