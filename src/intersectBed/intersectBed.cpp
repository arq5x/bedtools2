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
        // -wao the user wants to force the reporting of 0 overlap
        if (hits.size() == 0) {
            if (_writeAllOverlap) { 
                _bedA->reportBedTab(a);
                _bedB->reportNullBedTab();
                printf("0\n");
            }
            else if (_leftJoin) {
                _bedA->reportBedTab(a);
                _bedB->reportNullBedNewLine();
            }
        }
        else {
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
                           float overlapFraction, bool noHit, bool leftJoin, bool writeCount, bool sameStrand, bool diffStrand,
                           bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput, bool isUncompressedBam,
                           bool sortedInput, bool printHeader) {

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _anyHit              = anyHit;
    _noHit               = noHit;
    _leftJoin            = leftJoin;
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


bool BedIntersect::FindBlockedOverlaps(const BED &a, const vector<BED> &a_blocks, 
                                       const vector<BED> &hits, bool a_is_bam) {
    int a_footprint = GetTotalBlockLength(a_blocks);
    // container to store the set of raw hits 
    // that actually overlap the A blocks
    bedVector valid_hits;
    valid_hits.reserve(hits.size());
    
    // 1. Loop through each raw hit (outer loop)
    // 2. Break the raw hit into it;s blocks
    //    and see of one of the hit blocks (inner loop)
    //    overlaps one of a's blocks (inner, inner loop)
    // 3. If so, mark the hit as valid and add it to the valid_set.
    //    Otherwise, the hit only overlapped the span of a and not
    //    the individual blocks.  Thus, it doesn't count.
    bedVector::const_iterator hItr = hits.begin();
    bedVector::const_iterator hEnd = hits.end();
    for (; hItr != hEnd; ++hItr) {
        // break the hit into blocks
        bedVector hitBlocks;
        GetBedBlocks(*hItr, hitBlocks);
        int b_footprint = GetTotalBlockLength(hitBlocks);
        // test to see if there is a valid hit with one of the blocks
        bool valid_hit    = false;
        int total_overlap = 0;
        bedVector::const_iterator hbItr = hitBlocks.begin();
        bedVector::const_iterator hbEnd = hitBlocks.end();
        for (; hbItr != hbEnd; ++hbItr) {
            // look for overlaps between this hit/block and each block in a
            bedVector::const_iterator a_blockItr = a_blocks.begin();
            bedVector::const_iterator a_blockEnd = a_blocks.end();
            for (; a_blockItr != a_blockEnd; ++a_blockItr) {
                int hs = max(a_blockItr->start, hbItr->start);
                int he = min(a_blockItr->end, hbItr->end);
                int overlap = he - hs;
                if (overlap > 0) {
                    valid_hit = true;
                    total_overlap += overlap;
                }
            }
        }
        if (valid_hit) {
            // require sufficint overlap fraction (reciprocal or otherwise)
            // w.r.t to the "footprint" (i.e., the total length of each block)
            if ( ((float) total_overlap / (float) a_footprint) > _overlapFraction) {
                if (_reciprocal && ((float) total_overlap / (float) b_footprint) > _overlapFraction) {
                    valid_hits.push_back(*hItr);
                }
                else if (!_reciprocal) {
                    valid_hits.push_back(*hItr);
                }
            }
        }
    }
    if (!a_is_bam) {
        return processHits(a, valid_hits);
    }
    else
        return !valid_hits.empty();
}


void BedIntersect::ReportOverlapDetail(int overlapBases, const BED &a, const BED &b, CHRPOS s, CHRPOS e) {
    // default. simple intersection only
    if (_writeA == false && _writeB == false && 
        _writeOverlap == false && _leftJoin == false) 
    {
        _bedA->reportBedRangeNewLine(a,s,e);
    }
    //  -wa -wb write the original A and B
    else if ((_writeA == true && _writeB == true) || _leftJoin == true) {
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
                if (_obeySplits == false)
                    FindOverlaps(a, hits);
                // split the BED12 into blocks and look for overlaps in each discrete block
                else {
                    // find the hits that overlap with the full span of the blocked BED
                    _bedB->allHits(a.chrom, a.start, a.end, a.strand,
                                   hits, _sameStrand, _diffStrand,
                                   _overlapFraction, _reciprocal);
                    // break a into discrete blocks, as we need to 
                    // measure overlap with the individual blocks, not the full span.
                    bedVector a_blocks; 
                    GetBedBlocks(a, a_blocks);
                    // find the overlaps between the block in A and B 
                    // last parm is false as a is not a BAM entry
                    FindBlockedOverlaps(a, a_blocks, hits, false);
                }
                hits.clear();
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
            if (_obeySplits == false)
                processHits(hit_set.first, hit_set.second);
            else {
                bedVector a_blocks; 
                GetBedBlocks(hit_set.first, a_blocks);
                FindBlockedOverlaps(hit_set.first, a_blocks, hit_set.second, false);
            }
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
        _bedA->bedType = 12;
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

        // save an unaligned read if -v
        if (!bam.IsMapped()) {
            if (_noHit == true)
                writer.SaveAlignment(bam);
            continue;
        }   
        // break alignment into discrete blocks,
        bedVector bed_blocks;
        string chrom = refs.at(bam.RefID).RefName;
        GetBamBlocks(bam, chrom, bed_blocks, false, true);
        // create a basic BED entry from the BAM alignment
        BED bed;
        MakeBedFromBam(bam, chrom, bed_blocks, bed);
        bool overlapsFound = false;
        if ((_bamOutput == true) && (_obeySplits == false))
        {
            overlapsFound = _bedB->anyHits(bed.chrom, bed.start, bed.end, 
                                           bed.strand, _sameStrand, _diffStrand,
                                           _overlapFraction, _reciprocal);
        }
        else if ( ((_bamOutput == true)  && (_obeySplits == true)) ||
                  ((_bamOutput == false) && (_obeySplits == true)) )
        {
            // find the hits that overlap with the full span of the blocked BED
            _bedB->allHits(bed.chrom, bed.start, bed.end, bed.strand,
                           hits, _sameStrand, _diffStrand,
                           _overlapFraction, _reciprocal);
            // find the overlaps between the block in A and B
            overlapsFound = FindBlockedOverlaps(bed, bed_blocks, hits, _bamOutput);
        }
        else if ((_bamOutput == false) && (_obeySplits == false))
        {
            FindOverlaps(bed, hits);
        }
        // save the BAM alignment if overlap reqs. were met
        if (_bamOutput == true) {
            if ((overlapsFound == true) && (_noHit == false))
                writer.SaveAlignment(bam);
            else if ((overlapsFound == false) && (_noHit == true))
                writer.SaveAlignment(bam);
        }
        hits.clear();
    }

    // close the relevant BAM files.
    reader.Close();
    if (_bamOutput == true) {
        writer.Close();
    }
}

