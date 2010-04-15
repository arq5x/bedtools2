/*****************************************************************************
  intersectBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "intersectBed.h"


/*
	Constructor
*/
BedIntersect::BedIntersect(string bedAFile, string bedBFile, bool anyHit, 
						   bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
						   float overlapFraction, bool noHit, bool writeCount, bool forceStrand, 
						   bool reciprocal, bool bamInput, bool bamOutput) {

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
	_bamInput            = bamInput;
	_bamOutput           = bamOutput;
	
	// create new BED file objects for A and B
	_bedA = new BedFile(bedAFile);
	_bedB = new BedFile(bedBFile);
}


/*
	Destructor
*/
BedIntersect::~BedIntersect(void) {
}


bool BedIntersect::FindOverlaps(const BED &a, vector<BED> &hits) {
	
	bool hitsFound = false;
	
	// grab _all_ of the features in B that overlap with a.
	_bedB->FindOverlapsPerBin(a.chrom, a.start, a.end, a.strand, hits, _forceStrand); 
	
	// how many overlaps are there b/w a and B?
	int numOverlaps = 0;		
	
	// should we print each overlap, or does the user want summary information?
	bool printable = true;			
	if (_anyHit || _noHit || _writeCount)
		printable = false;
	
	// loop through the hits and report those that meet the user's criteria
	vector<BED>::const_iterator h = hits.begin();
	vector<BED>::const_iterator hitsEnd = hits.end();
	for (; h != hitsEnd; ++h) {
		
		int s            = max(a.start, h->start);
		int e            = min(a.end, h->end);
		int overlapBases = (e - s);				// the number of overlapping bases b/w a and b
		int aLength      = (a.end - a.start);   // the length of a in b.p.
		
		// is there enough overlap relative to the user's request? (default ~ 1bp)
		if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) { 
		
			// Report the hit if the user doesn't care about reciprocal overlap between A and B.
			if (_reciprocal == false) {
				hitsFound = true;
				numOverlaps++;
			}
			// we require there to be sufficient __reciprocal__ overlap
			else {			
				int bLength    = (h->end - h->start);
				float bOverlap = ( (float) overlapBases / (float) bLength );
				if (bOverlap >= _overlapFraction) {
					hitsFound = true;
					numOverlaps++;
				}
			}
			// report the overlap per the user's request.
			if (printable == true) {
				ReportOverlapDetail(overlapBases, a, *h, s, e);
			}
		}
	}
	// report the summary of the overlaps if requested.
	ReportOverlapSummary(a, numOverlaps);
	// were hits found for this BED feature?
	return hitsFound;
}


void BedIntersect::ReportOverlapDetail(const int &overlapBases, const BED &a, const BED &b,
									   const int &s, const int &e) {
	// simple intersection only
	if (_writeA == false && _writeB == false && _writeOverlap == false) {
		_bedA->reportBedRangeNewLine(a,s,e);
	}
	// write the original A and B 
	else if (_writeA == true && _writeB == true) {
		_bedA->reportBedTab(a);
		_bedB->reportBedNewLine(b);
	}
	// write just the original A
	else if (_writeA == true) {
		_bedA->reportBedNewLine(a);
	}
	// write the intersected portion of A and the original B 
	else if (_writeB == true) {
		_bedA->reportBedRangeTab(a,s,e);
		_bedB->reportBedNewLine(b);
	}
	// write the original A and B plus the no. of overlapping bases.
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
	BED a;	
	
	// open the "A" file, process each BED entry and searh for overlaps.
	_bedA->Open();
	while (_bedA->GetNextBed(a, lineNum) == true) {
		FindOverlaps(a, hits);
		hits.clear();
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
		writer.Open("stdout", header, refs);
	}

	vector<BED> hits;					// vector of potential hits
	// reserve some space
	hits.reserve(100);
	
	_bedA->bedType = 6;
	BamAlignment bam;	
	bool overlapsFound;
	// get each set of alignments for each pair.
	while (reader.GetNextAlignment(bam)) {
		
		if (bam.IsMapped()) {	
			BED a;
			a.chrom = refs.at(bam.RefID).RefName;
			a.start = bam.Position;
			a.end   = bam.Position + bam.AlignedBases.size();

			// build the name field from the BAM alignment.
			a.name = bam.Name;
			if (bam.IsFirstMate()) a.name += "/1";
			if (bam.IsSecondMate()) a.name += "/2";

			a.score  = ToString(bam.MapQuality);
			a.strand = "+"; if (bam.IsReverseStrand()) a.strand = "-"; 
	
			if (_bamOutput == true) {
				overlapsFound = FindOneOrMoreOverlap(a);
				if (overlapsFound == true) {
					if (_noHit == false) 
						writer.SaveAlignment(bam);
				}
				else {
					if (_noHit == true) 
						writer.SaveAlignment(bam);	
				}
			}
			else {
				overlapsFound = FindOverlaps(a, hits);
				hits.clear();
			}
		}
	}
	
	// close the relevant BAM files.
	reader.Close();
	if (_bamOutput == true) {
		writer.Close();
	}
}


void BedIntersect::DetermineBedInput() {

	// dealing with a proper file
	if (_bedA->bedFile != "stdin") {   
		// it's BED or GFF
		if (_bamInput == false) { 
			IntersectBed();
		}
		else {
			IntersectBam(_bedA->bedFile);
		}
	}
	// reading from stdin
	else {  
		// it's BED or GFF 
		if (_bamInput == false) {					
			IntersectBed();
		}
		// it's BAM
		else {
			IntersectBam("stdin");
		}				
	}
}
