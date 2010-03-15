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
						   bool writeA, bool writeB, bool writeOverlap, float overlapFraction, 
						   bool noHit, bool writeCount, bool forceStrand, bool reciprocal,
						   bool bamInput, bool bamOutput) {

	this->bedAFile        = bedAFile;
	this->bedBFile        = bedBFile;
	this->anyHit          = anyHit;
	this->noHit           = noHit;
	this->writeA          = writeA;	
	this->writeB          = writeB;
	this->writeOverlap    = writeOverlap;	
	this->writeCount      = writeCount;
	this->overlapFraction = overlapFraction;
	this->forceStrand     = forceStrand;
	this->reciprocal      = reciprocal;
	this->bamInput        = bamInput;
	this->bamOutput       = bamOutput;
	
	// create new BED file objects for A and B
	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
}


/*
	Destructor
*/
BedIntersect::~BedIntersect(void) {
}


bool BedIntersect::FindOverlaps(const BED &a, vector<BED> &hits) {
	
	bool hitsFound = false;
	
	// grab _all_ of the features in B that overlap with a.
	bedB->FindOverlapsPerBin(a.chrom, a.start, a.end, a.strand, hits, this->forceStrand); 
	
	// how many overlaps are there b/w a and B?
	int numOverlaps = 0;		
	
	// should we print each overlap, or does the user want summary information?
	bool printable = true;			
	if (anyHit || noHit || writeCount) {
		printable = false;
	}
	
	// loop through the hits and report those that meet the user's criteria
	vector<BED>::const_iterator h = hits.begin();
	vector<BED>::const_iterator hitsEnd = hits.end();
	for (; h != hitsEnd; ++h) {
		
		int s = max(a.start, h->start);
		int e = min(a.end, h->end);
		int overlapBases = (e - s);				// the number of overlapping bases b/w a and b
		int aLength = (a.end - a.start);		// the length of a in b.p.
		
		// is there enough overlap relative to the user's request? (default ~ 1bp)
		if ( ( (float) overlapBases / (float) aLength ) >= this->overlapFraction ) { 
		
			// Report the hit if the user doesn't care about reciprocal overlap between A and B.
			if (!reciprocal) {
				
				hitsFound = true;
				numOverlaps++;		// we have another hit for A
/*				
				if (!writeB && printable) {
					if (writeA) bedA->reportBedNewLine(a); 
					else bedA->reportBedRangeNewLine(a,s,e);
				}
				else if (printable) {
					if (writeA) {
						bedA->reportBedTab(a);
						bedB->reportBedNewLine(*h);
					}
					else {
						bedA->reportBedRangeTab(a,s,e);
						bedB->reportBedNewLine(*h);									
					}
				}
*/				
				// new
				if (printable == true) {
					if (writeA == false && writeB == false && writeOverlap == false) {
						bedA->reportBedRangeNewLine(a,s,e);
					}
					else if (writeA == true && writeB == true) {
						bedA->reportBedTab(a);
						bedB->reportBedNewLine(*h);
					}
					else if (writeA == true) {
						bedA->reportBedNewLine(a);
					}
					else if (writeB == true) {
						bedA->reportBedRangeTab(a,s,e);
						bedB->reportBedNewLine(*h);
					}
					else if (writeOverlap == true) {
						bedA->reportBedTab(a);
						bedB->reportBedTab(*h);
						printf("%d\n", overlapBases);
					}
				}
				
			}
			else {			// the user wants there to be sufficient reciprocal overlap
				int bLength = (h->end - h->start);
				float bOverlap = ( (float) overlapBases / (float) bLength );
			
				if (bOverlap >= this->overlapFraction) {
					hitsFound = true;
					numOverlaps++;		// we have another hit for A
				
					if (!writeB && printable) {
						if (writeA) bedA->reportBedNewLine(a);
						else bedA->reportBedRangeNewLine(a,s,e);
					}
					else if (printable) {
						if (writeA) {
							bedA->reportBedTab(a);
							bedB->reportBedNewLine(*h);
						}
						else {
							bedA->reportBedRangeTab(a,s,e);
							bedB->reportBedNewLine(*h);									
						}
					}
				}
			}
		}
	}
	if (anyHit && (numOverlaps >= 1)) {
		bedA->reportBedNewLine(a);
	}
	else if (writeCount) {
		bedA->reportBedTab(a); 
		printf("%d\n", numOverlaps);
	}
	else if (noHit && (numOverlaps == 0)) {
		bedA->reportBedNewLine(a);
	}

	return hitsFound;
}


bool BedIntersect::FindOneOrMoreOverlap(const BED &a) {
	bool overlapsFound;
	if (this->reciprocal == false) {
		overlapsFound = bedB->FindOneOrMoreOverlapsPerBin(a.chrom, a.start, a.end, a.strand, 
			this->forceStrand, this->overlapFraction); 
	}
	else {
		overlapsFound = bedB->FindOneOrMoreReciprocalOverlapsPerBin(a.chrom, a.start, a.end, a.strand, 
			this->forceStrand, this->overlapFraction);
	}
	return overlapsFound;
}
 

void BedIntersect::IntersectBed(istream &bedInput) {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();
	string bedLine;                                                                                                                    
	int lineNum = 0;					// current input line number
	vector<BED> hits;					// vector of potential hits
	vector<string> bedFields;			// vector for a BED entry
	
	// reserve some space
	hits.reserve(100);
	bedFields.reserve(12);	
		
	// process each entry in A
	while (getline(bedInput, bedLine)) {

		lineNum++;
		Tokenize(bedLine,bedFields);
		BED a;
		// find the overlaps with B if it's a valid BED entry. 
		if (bedA->parseLine(a, bedFields, lineNum)) {
			FindOverlaps(a, hits);
			hits.clear();
		}
		// reset for the next input line
		bedFields.clear();
	}
}


void BedIntersect::IntersectBam(string bamFile) {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();
	
	// open the BAM file
	BamReader reader;
	BamWriter writer;
	reader.Open(bamFile);

	// get header & reference information
	string header = reader.GetHeaderText();
	RefVector refs = reader.GetReferenceData();

	// open a BAM output to stdout if we are writing BAM
	if (this->bamOutput == true) {
		// open our BAM writer
		writer.Open("stdout", header, refs);
	}

	vector<BED> hits;					// vector of potential hits
	// reserve some space
	hits.reserve(100);
	
	bedA->bedType = 6;
	BamAlignment bam;	
	bool overlapsFound;
	// get each set of alignments for each pair.
	while (reader.GetNextAlignment(bam)) {
		
		if (bam.IsMapped()) {	
			BED a;
			a.chrom = refs.at(bam.RefID).RefName;
			a.start = bam.Position;
			a.end = bam.Position + bam.AlignedBases.size();

			// build the name field from the BAM alignment.
			a.name = bam.Name;
			if (bam.IsFirstMate()) a.name += "/1";
			if (bam.IsSecondMate()) a.name += "/2";

			a.score = ToString(bam.MapQuality);
			a.strand = "+"; if (bam.IsReverseStrand()) a.strand = "-"; 
	
			if (this->bamOutput == true) {
				overlapsFound = FindOneOrMoreOverlap(a);
				if (overlapsFound == true) {
					if (!this->noHit) writer.SaveAlignment(bam);
				}
				else {
					if (this->noHit) writer.SaveAlignment(bam);	
				}
			}
			else {
				overlapsFound = FindOverlaps(a, hits);
				hits.clear();
			}
		}
	}
	
	reader.Close();
	if (this->bamOutput == true) {
		writer.Close();
	}
}


void BedIntersect::DetermineBedInput() {
	
	if (bedA->bedFile != "stdin") {   // process a file
		if (this->bamInput == false) { //bed/gff
			ifstream beds(bedA->bedFile.c_str(), ios::in);
			if ( !beds ) {
				cerr << "Error: The requested bed file (" << bedA->bedFile << ") could not be opened. Exiting!" << endl;
				exit (1);
			}
			IntersectBed(beds);
		}
		else {	// bam
			IntersectBam(bedA->bedFile);
		}
	}
	else {   // process stdin
		if (this->bamInput == false) {	//bed/gff					
			IntersectBed(cin);
		}
		else {
			IntersectBam("stdin");
		}				
	}
}
