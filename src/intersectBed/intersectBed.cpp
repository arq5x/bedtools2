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
bool BedIntersect::processHits(const BED &a, const vector<BED> &hits) {

    bool hitsFound   = false;
    if (_printable == true) {
        vector<BED>::const_iterator h       = hits.begin();
        vector<BED>::const_iterator hitsEnd = hits.end();
        for (; h != hitsEnd; ++h) {
            CHRPOS s = max(a.start, h->start);
            CHRPOS e = min(a.end, h->end);
            int overlapBases = (e - s);
            ReportOverlapDetail(overlapBases, a, *h, s, e);
            hitsFound = true;
        }
    }
    else {
        ReportOverlapSummary(a, hits.size());
    }
    return hitsFound;
}


/*
    Constructor
*/
BedIntersect::BedIntersect(string bedAFile, string bedBFile, bool anyHit,
                           bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                           float overlapFraction, bool noHit, bool writeCount, bool sameStrand, bool diffStrand,
                           bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput, bool isUncompressedBam,
                           bool sortedInput, bool printHeader) {

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
    _sameStrand          = sameStrand;
    _diffStrand          = diffStrand;
    _reciprocal          = reciprocal;
    _obeySplits          = obeySplits;
    _bamInput            = bamInput;
    _bamOutput           = bamOutput;
    _isUncompressedBam   = isUncompressedBam;
    _sortedInput         = sortedInput;
    _printHeader         = printHeader;

    // should we print each overlap, or does the user want summary information?
    _printable = true;
    if (_anyHit || _noHit || _writeCount)
        _printable = false;
        
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
    
    // collect and report the sufficient hits
    _bedB->allHits(a.chrom, a.start, a.end, a.strand,
                   hits, _sameStrand, _diffStrand,
                   _overlapFraction, _reciprocal);
    hitsFound = processHits(a, hits);
    return hitsFound;
}


void BedIntersect::ReportOverlapDetail(int overlapBases, const BED &a, const BED &b, CHRPOS s, CHRPOS e) {
    // default. simple intersection only
    if (_writeA == false && _writeB == false && _writeOverlap == false) {
        _bedA->reportBedRangeNewLine(a,s,e);
    }
    //  -wa -wb write the original A and B
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
        if (b.zeroLength == false) printf("%d\n", overlapBases);
        else printf("0\n");
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


void BedIntersect::IntersectBed() {

    // create new BED file objects for A and B
    _bedA = new BedFile(_bedAFile);
    _bedB = new BedFile(_bedBFile);

    if (_sortedInput == false) {
        // load the "B" file into a map in order to
        // compare each entry in A to it in search of overlaps.
        _bedB->loadBedFileIntoMap();

        vector<BED> hits;
        hits.reserve(100);
        BED a;

        // open the "A" file, process each BED entry and searh for overlaps.
        _bedA->Open();
        // report A's header first if asked.
        if (_printHeader == true) {
            _bedA->PrintHeader();
        }
        while (_bedA->GetNextBed(a)) {
            if (_bedA->_status == BED_VALID) {
                // treat the BED as a single "block"
                if (_obeySplits == false) {
                    FindOverlaps(a, hits);
                    hits.clear();
                }
                // split the BED12 into blocks and look for overlaps in each discrete block
                else {
                    bedVector bedBlocks;  // vec to store the discrete BED "blocks"
                    GetBedBlocks(a,bedBlocks);

                    vector<BED>::const_iterator bedItr  = bedBlocks.begin();
                    vector<BED>::const_iterator bedEnd  = bedBlocks.end();
                    for (; bedItr != bedEnd; ++bedItr) {
                        FindOverlaps(*bedItr, hits);
                        hits.clear();
                    }
                }
            }
        }
        _bedA->Close();
    }
    else {
        // use the chromsweep algorithm to detect overlaps on the fly.
        ChromSweep sweep = ChromSweep(_bedA, _bedB, 
                                      _sameStrand, _diffStrand, 
                                      _overlapFraction, _reciprocal,
                                      _printHeader);

        pair<BED, vector<BED> > hit_set;
        hit_set.second.reserve(10000);
        while (sweep.Next(hit_set)) {
            processHits(hit_set.first, hit_set.second);
        }
    }
}


void BedIntersect::IntersectBam(string bamFile) {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB = new BedFile(_bedBFile);
    _bedB->loadBedFileIntoMap();

    // create a dummy BED A file for printing purposes if not
    // using BAM output.
    if (_bamOutput == false) {
        _bedA = new BedFile(_bedAFile);
        _bedA->bedType = 6;
    }

    // open the BAM file
    BamReader reader;
    BamWriter writer;
    
    reader.Open(bamFile);

    // get header & reference information
    string bamHeader  = reader.GetHeaderText();
    RefVector refs    = reader.GetReferenceData();

    // open a BAM output to stdout if we are writing BAM
    if (_bamOutput == true) {
        // set compression mode
        BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
        if ( _isUncompressedBam ) compressionMode = BamWriter::Uncompressed;
        writer.SetCompressionMode(compressionMode);
        // open our BAM writer
        writer.Open("stdout", bamHeader, refs);
    }

    vector<BED> hits;
    // reserve some space
    hits.reserve(100);

    
    BamAlignment bam;    
    // get each set of alignments for each pair.
    while (reader.GetNextAlignment(bam)) {

        if (bam.IsMapped()) {
            BED a;
            a.chrom = refs.at(bam.RefID).RefName;
            a.start = bam.Position;
            a.end   = bam.GetEndPosition(false, false);

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
                    overlapsFound = _bedB->anyHits(a.chrom, a.start, a.end, 
                                                   a.strand, _sameStrand, _diffStrand,
                                                   _overlapFraction, _reciprocal);
                }
                // split the BAM alignment into discrete blocks and
                // look for overlaps only within each block.
                else {
                    bool overlapFoundForBlock;
                    bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
                    // we don't want to split on "D" ops, hence the "false"
                    GetBamBlocks(bam, refs.at(bam.RefID).RefName, bedBlocks);

                    vector<BED>::const_iterator bedItr  = bedBlocks.begin();
                    vector<BED>::const_iterator bedEnd  = bedBlocks.end();
                    for (; bedItr != bedEnd; ++bedItr) {
                        overlapFoundForBlock = _bedB->anyHits(bedItr->chrom, bedItr->start, bedItr->end, 
                                                              bedItr->strand, _sameStrand, _diffStrand,
                                                              _overlapFraction, _reciprocal);
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
                    GetBamBlocks(bam, refs.at(bam.RefID).RefName, bedBlocks);

                    vector<BED>::const_iterator bedItr  = bedBlocks.begin();
                    vector<BED>::const_iterator bedEnd  = bedBlocks.end();
                    for (; bedItr != bedEnd; ++bedItr) {
                        FindOverlaps(*bedItr, hits);
                        hits.clear();
                    }
                }
            }
        }
        // BAM IsMapped() is false
        else if (_noHit == true) {
            writer.SaveAlignment(bam);
        }
    }

    // close the relevant BAM files.
    reader.Close();
    if (_bamOutput == true) {
        writer.Close();
    }
}

