/*****************************************************************************
  pairToBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "pairToBed.h"


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


/*
    Constructor
*/


BedIntersectPE::BedIntersectPE(string bedAFilePE, string bedBFile, float overlapFraction,
                               string searchType, bool sameStrand, bool diffStrand, bool bamInput,
                               bool bamOutput,  bool uncompressedBam, bool useEditDistance) {

    _bedAFilePE        = bedAFilePE;
    _bedBFile          = bedBFile;
    _overlapFraction   = overlapFraction;
    _sameStrand        = sameStrand;
    _diffStrand        = diffStrand;
    _useEditDistance   = useEditDistance;
    _searchType        = searchType;
    _bamInput          = bamInput;
    _bamOutput         = bamOutput;
    _isUncompressedBam = uncompressedBam;

    _bedA = new BedFilePE(bedAFilePE);
    _bedB = new BedFile(bedBFile);

    if (_bamInput == false)
        IntersectBedPE();
    else
        IntersectBamPE(_bedAFilePE);
}


/*
    Destructor
*/

BedIntersectPE::~BedIntersectPE(void) {
}



void BedIntersectPE::FindOverlaps(const BEDPE &a, vector<BED> &hits1, vector<BED> &hits2, const string &type) {

    // list of hits on each end of BEDPE
    // that exceed the requested overlap fraction
    vector<BED> qualityHits1;
    vector<BED> qualityHits2;

    // count of hits on each end of BEDPE
    // that exceed the requested overlap fraction
    int numOverlapsEnd1 = 0;
    int numOverlapsEnd2 = 0;

    // make sure we have a valid chromosome before we search
    if (a.chrom1 != ".") {
        // Find the quality hits between ***end1*** of the BEDPE and the B BED file
        _bedB->allHits(a.chrom1, a.start1, a.end1, a.strand1, 
                       hits1, _sameStrand, _diffStrand, 0.0, false);

        vector<BED>::const_iterator h = hits1.begin();
        vector<BED>::const_iterator hitsEnd = hits1.end();
        for (; h != hitsEnd; ++h) {

            int s = max(a.start1, h->start);
            int e = min(a.end1, h->end);
            int overlapBases = (e - s);             // the number of overlapping bases b/w a and b
            int aLength = (a.end1 - a.start1);      // the length of a in b.p.

            // is there enough overlap relative to the user's request? (default ~ 1bp)
            if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) {
                numOverlapsEnd1++;

                if (type == "either") {
                    _bedA->reportBedPETab(a);
                    _bedB->reportBedNewLine(*h);
                }
                else {
                    qualityHits1.push_back(*h);
                }
            }
        }
    }


    // make sure we have a valid chromosome before we search
    if (a.chrom2 != ".") {
        // Now find the quality hits between ***end2*** of the BEDPE and the B BED file
        _bedB->allHits(a.chrom2, a.start2, a.end2, a.strand2, 
                       hits2, _sameStrand, _diffStrand, 0.0, false);

        vector<BED>::const_iterator h = hits2.begin();
        vector<BED>::const_iterator hitsEnd = hits2.end();
        for (; h != hitsEnd; ++h) {

            int s = max(a.start2, h->start);
            int e = min(a.end2, h->end);
            int overlapBases = (e - s);             // the number of overlapping bases b/w a and b
            int aLength = (a.end2 - a.start2);      // the length of a in b.p.

            // is there enough overlap relative to the user's request? (default ~ 1bp)
            if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) {
                numOverlapsEnd2++;

                if (type == "either") {
                    _bedA->reportBedPETab(a);
                    _bedB->reportBedNewLine(*h);
                }
                else {
                    qualityHits2.push_back(*h);
                }
            }
        }
    }

    // Now report the hits depending on what the user has requested.
    if (type == "neither") {
        if ( (numOverlapsEnd1 == 0) && (numOverlapsEnd2 == 0) ) {
            _bedA->reportBedPENewLine(a);
        }
    }
    else if (type == "notboth") {
        if ( (numOverlapsEnd1 == 0) && (numOverlapsEnd2 == 0) ) {
            _bedA->reportBedPENewLine(a);
        }
        else if ( (numOverlapsEnd1 > 0) && (numOverlapsEnd2 == 0) ) {
            for (vector<BED>::iterator q = qualityHits1.begin(); q != qualityHits1.end(); ++q) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*q);
            }
        }
        else if ( (numOverlapsEnd1 == 0) && (numOverlapsEnd2 > 0) ) {
            for (vector<BED>::iterator q = qualityHits2.begin(); q != qualityHits2.end(); ++q) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*q);
            }
        }
    }
    else if (type == "xor") {
        if ( (numOverlapsEnd1 > 0) && (numOverlapsEnd2 == 0) ) {
            for (vector<BED>::iterator q = qualityHits1.begin(); q != qualityHits1.end(); ++q) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*q);
            }
        }
        else if ( (numOverlapsEnd1 == 0) && (numOverlapsEnd2 > 0) ) {
            for (vector<BED>::iterator q = qualityHits2.begin(); q != qualityHits2.end(); ++q) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*q);
            }
        }
    }
    else if (type == "both") {
        if ( (numOverlapsEnd1 > 0) && (numOverlapsEnd2 > 0) ) {
            for (vector<BED>::iterator q1 = qualityHits1.begin(); q1 != qualityHits1.end(); ++q1) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*q1);
            }
            for (vector<BED>::iterator q2 = qualityHits2.begin(); q2 != qualityHits2.end(); ++q2) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*q2);
            }
        }
    }
}


bool BedIntersectPE::FindOneOrMoreOverlaps(const BEDPE &a, const string &type) {

    // flags for the existence of hits on each end of BEDPE
    // that exceed the requested overlap fraction
    bool end1Found = false;
    bool end2Found = false;

    // Look for overlaps in end 1 assuming we have an aligned chromosome.
    if (a.chrom1 != ".") {
        end1Found = _bedB->anyHits(a.chrom1, a.start1, a.end1, a.strand1,
                                   _sameStrand, _diffStrand, _overlapFraction, false);

        // can we bail out without checking end2?
        if ((type == "either") && (end1Found == true)) return true;
        else if ((type == "neither") && (end1Found == true)) return false;
        else if ((type == "notboth") && (end1Found == false)) return true;
        else if ((type == "both") && (end1Found == false)) return false;
    }

    // Now look for overlaps in end 2 assuming we have an aligned chromosome.
    if (a.chrom2 != ".") {
        end2Found = _bedB->anyHits(a.chrom2, a.start2, a.end2, a.strand2,
                                   _sameStrand, _diffStrand, _overlapFraction, false);

        if ((type == "either") && (end2Found == true)) return true;
        else if ((type == "neither") && (end2Found == true)) return false;
        else if ((type == "notboth") && (end2Found == false)) return true;
        else if ((type == "both") && (end2Found == false)) return false;
    }

    // Now report the hits depending on what the user has requested.
    if (type == "notboth") {
        if ( (end1Found == false) || (end2Found == false) ) return true;
        else return false;
    }
    else if (type == "either") {
        if ( (end1Found == false) && (end2Found == false) ) return false;
    }
    else if (type == "neither") {
        if ( (end1Found == false) && (end2Found == false) ) return true;
        else return false;
    }
    else if (type == "xor") {
        if ( (end1Found == true) && (end2Found == false) ) return true;
        else if ( (end1Found == false) && (end2Found == true) ) return true;
        else return false;
    }
    else if (type == "both") {
        if ( (end1Found == true) && (end2Found == true) ) return true;
        return false;
    }
    return false;
}


void BedIntersectPE::FindSpanningOverlaps(const BEDPE &a, vector<BED> &hits, const string &type) {

    // count of hits on _between_ end of BEDPE
    // that exceed the requested overlap fraction
    int numOverlaps = 0;
    CHRPOS spanStart = 0;
    CHRPOS spanEnd = 0;
    CHRPOS spanLength = 0;

    if ((type == "ispan") || (type == "notispan")) {
        spanStart = a.end1;
        spanEnd = a.start2;
        if (a.end1 > a.start2) {
            spanStart = a.end2;
            spanEnd = a.start1;
        }
    }
    else if ((type == "ospan") || (type == "notospan")) {
        spanStart = a.start1;
        spanEnd = a.end2;
        if (a.start1 > a.start2) {
            spanStart = a.start2;
            spanEnd = a.end1;
        }
    }
    spanLength = spanEnd - spanStart;

    // get the hits for the span
    _bedB->allHits(a.chrom1, spanStart, spanEnd, a.strand1, 
                   hits, _sameStrand, _diffStrand, 0.0, false);

    vector<BED>::const_iterator h = hits.begin();
    vector<BED>::const_iterator hitsEnd = hits.end();
    for (; h != hitsEnd; ++h) {

        CHRPOS s = max(spanStart, h->start);
        CHRPOS e = min(spanEnd, h->end);
        CHRPOS overlapBases = (e - s);                     // the number of overlapping bases b/w a and b
        CHRPOS spanLength = (spanEnd - spanStart);     // the length of a in b.p.

        // is there enough overlap relative to the user's request? (default ~ 1bp)
        if ( ( (float) overlapBases / (float) spanLength ) >= _overlapFraction ) {
            numOverlaps++;
            if ((type == "ispan") || (type == "ospan")) {
                _bedA->reportBedPETab(a);
                _bedB->reportBedNewLine(*h);
            }
        }
    }

    if ( ( (type == "notispan") || (type == "notospan") ) && numOverlaps == 0 ) {
        _bedA->reportBedPENewLine(a);
    }
}


bool BedIntersectPE::FindOneOrMoreSpanningOverlaps(const BEDPE &a, const string &type) {

    CHRPOS spanStart = 0;
    CHRPOS spanEnd = 0;
    CHRPOS spanLength = 0;
    bool overlapFound;

    if ((type == "ispan") || (type == "notispan")) {
        spanStart = a.end1;
        spanEnd = a.start2;
        if (a.end1 > a.start2) {
            spanStart = a.end2;
            spanEnd = a.start1;
        }
    }
    else if ((type == "ospan") || (type == "notospan")) {
        spanStart = a.start1;
        spanEnd = a.end2;
        if (a.start1 > a.start2) {
            spanStart = a.start2;
            spanEnd = a.end1;
        }
    }
    spanLength = spanEnd - spanStart;

    overlapFound = _bedB->anyHits(a.chrom1, spanStart, spanEnd, a.strand1,
                                  _sameStrand, _diffStrand, _overlapFraction, false);

    return overlapFound;
}


void BedIntersectPE::IntersectBedPE() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedFileIntoMap();

    int lineNum = 0;                    // current input line number
    vector<BED> hits, hits1, hits2;     // vector of potential hits

    // reserve some space
    hits.reserve(100);
    hits1.reserve(100);
    hits2.reserve(100);

    BEDPE a, nullBedPE;
    BedLineStatus bedStatus;

    _bedA->Open();
    while ((bedStatus = _bedA->GetNextBedPE(a, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            if ( (_searchType == "ispan") || (_searchType == "ospan") ||
                 (_searchType == "notispan") || (_searchType == "notospan") ) {
                if (a.chrom1 == a.chrom2) {
                    FindSpanningOverlaps(a, hits, _searchType);
                    hits.clear();
                }
            }
            else {
                FindOverlaps(a, hits1, hits2, _searchType);
                hits1.clear();
                hits2.clear();
            }
            a = nullBedPE;
        }
    }
    _bedA->Close();
}


void BedIntersectPE::IntersectBamPE(string bamFile) {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedFileIntoMap();

    // open the BAM file
    BamReader reader;
    BamWriter writer;
    if (!reader.Open(bamFile)) {
        cerr << "Failed to open BAM file " << bamFile << endl;
        exit(1);
    }

    // get header & reference information
    string bamHeader = reader.GetHeaderText();
    RefVector refs   = reader.GetReferenceData();

    // open a BAM output to stdout if we are writing BAM
    if (_bamOutput == true) {
        // set compression mode
        BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
        if ( _isUncompressedBam ) compressionMode = BamWriter::Uncompressed;
        writer.SetCompressionMode(compressionMode);
        // open our BAM writer
        writer.Open("stdout", bamHeader, refs);
    }

    // track the previous and current sequence
    // names so that we can identify blocks of
    // alignments for a given read ID.
    string prevName, currName;
    prevName = currName = "";

    vector<BamAlignment> alignments;        // vector of BAM alignments for a given ID in a BAM file.
    alignments.reserve(100);

    _bedA->bedType = 10;                    // it's a full BEDPE given it's BAM

    // rip through the BAM file and convert each mapped entry to BEDPE
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        reader.GetNextAlignment(bam2);        
        if (bam1.Name != bam2.Name) {
            while (bam1.Name != bam2.Name)
            {
                if (bam1.IsPaired()) 
                {
                    cerr << "*****WARNING: Query " << bam1.Name
                         << " is marked as paired, but it's mate does not occur"
                         << " next to it in your BAM file.  Skipping. " << endl;
                }
                bam1 = bam2;
                reader.GetNextAlignment(bam2);
            }
        }
        else if (bam1.IsPaired() && bam1.IsPaired()) {
            ProcessBamBlock(bam1, bam2, refs, writer);
        }
    }
    // close up
    reader.Close();
    if (_bamOutput == true) {
        writer.Close();
    }
}


void BedIntersectPE::ProcessBamBlock (const BamAlignment &bam1, const BamAlignment &bam2,
                                      const RefVector &refs, BamWriter &writer) {

    vector<BED> hits, hits1, hits2;         // vector of potential hits
    hits.reserve(1000);                     // reserve some space
    hits1.reserve(1000);
    hits2.reserve(1000);

    bool overlapsFound;                     // flag to indicate if overlaps were found

    if ( (_searchType == "either") || (_searchType == "xor") ||
              (_searchType == "both") || (_searchType == "notboth") ||
              (_searchType == "neither") ) {

        // create a new BEDPE feature from the BAM alignments.
        BEDPE a;
        ConvertBamToBedPE(bam1, bam2, refs, a);
        if (_bamOutput == true) {   // BAM output
            // write to BAM if correct hits found
            overlapsFound = FindOneOrMoreOverlaps(a, _searchType);
            if (overlapsFound == true) {
                writer.SaveAlignment(bam1);
                writer.SaveAlignment(bam2);
            }
        }
        else {  // BEDPE output
            FindOverlaps(a, hits1, hits2, _searchType);
            hits1.clear();
            hits2.clear();
        }
    }
    else if ( (_searchType == "ispan") || (_searchType == "ospan") ) {
        // only look for ispan and ospan when both ends are mapped.
        if (bam1.IsMapped() && bam2.IsMapped()) {
            // only do an inspan or outspan check if the alignment is intrachromosomal
            if (bam1.RefID == bam2.RefID) {
                // create a new BEDPE feature from the BAM alignments.
                BEDPE a;
                ConvertBamToBedPE(bam1, bam2, refs, a);
                if (_bamOutput == true) {   // BAM output
                    // look for overlaps, and write to BAM if >=1 were found
                    overlapsFound = FindOneOrMoreSpanningOverlaps(a, _searchType);
                    if (overlapsFound == true) {
                        writer.SaveAlignment(bam1);
                        writer.SaveAlignment(bam2);
                    }
                }
                else {  // BEDPE output
                    FindSpanningOverlaps(a, hits, _searchType);
                    hits.clear();
                }
            }
        }
    }
    else if ( (_searchType == "notispan") || (_searchType == "notospan") ) {
        // only look for notispan and notospan when both ends are mapped.
        if (bam1.IsMapped() && bam2.IsMapped()) {
            // only do an inspan or outspan check if the alignment is intrachromosomal
            if (bam1.RefID == bam2.RefID) {
                // create a new BEDPE feature from the BAM alignments.
                BEDPE a;
                ConvertBamToBedPE(bam1, bam2, refs, a);
                if (_bamOutput == true) {   // BAM output
                    // write to BAM if there were no overlaps
                    overlapsFound = FindOneOrMoreSpanningOverlaps(a, _searchType);
                    if (overlapsFound == false) {
                        writer.SaveAlignment(bam1);
                        writer.SaveAlignment(bam2);
                    }
                }
                else {  // BEDPE output
                    FindSpanningOverlaps(a, hits, _searchType);
                    hits.clear();
                }
            }
            // if inter-chromosomal or orphaned, we know it's not ispan and not ospan
            else if (_bamOutput == true) {
                writer.SaveAlignment(bam1);
                writer.SaveAlignment(bam2);
            }
        }
        // if both ends aren't mapped, we know that it's notispan and not ospan
        else if (_bamOutput == true) {
            writer.SaveAlignment(bam1);
            writer.SaveAlignment(bam2);
        }
    }
}


