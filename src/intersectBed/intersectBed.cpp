#include "intersectBed.h"


BedIntersect::BedIntersect(string &bedAFile, string &bedBFile, bool &anyHit, 
bool &writeB, float &overlapFraction, bool &noHit, bool &writeCount) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	this->anyHit = anyHit;
	this->noHit = noHit;
	this->writeB = writeB;
	this->writeCount = writeCount;
	this->overlapFraction = overlapFraction;

	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
}



BedIntersect::~BedIntersect(void) {
}


void BedIntersect::reportA(const BED &a) {
	
	if (bedA->bedType == 3) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end;
	}
	else if (bedA->bedType == 4) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t"
		<< a.name;
	}
	else if (bedA->bedType == 5) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t"
		<< a.name << "\t" << a.score;
	}
	else if (bedA->bedType == 6) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t" 
		<< a.name << "\t" << a.score << "\t" << a.strand;
	}
	
}



void BedIntersect::reportAIntersect(const BED &a, int &start, int &end) {

	if (bedA->bedType == 3) {
		cout << a.chrom << "\t" << start << "\t" << end;
	}
	else if (bedA->bedType == 4) {
		cout << a.chrom << "\t" << start << "\t" << end << "\t"
		<< a.name;
	}
	else if (bedA->bedType == 5) {
		cout << a.chrom << "\t" << start << "\t" << end << "\t"
		<< a.name << "\t" << a.score;
	}
	else if (bedA->bedType == 6) {
		cout << a.chrom << "\t" << start << "\t" << end << "\t" 
		<< a.name << "\t" << a.score << "\t" << a.strand;
	}
	
}

void BedIntersect::reportB(const BED &b) {
	if (bedB->bedType == 3) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end;
	}
	else if (bedB->bedType == 4) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end << "\t"
		<< b.name;
	}
	else if (bedB->bedType == 5) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end << "\t"
		<< b.name << "\t" << b.score;
	}
	else if (bedB->bedType == 6) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end << "\t" 
		<< b.name << "\t" << b.score << "\t" << b.strand;
	}
}




 
void BedIntersect::IntersectBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();


	string bedLine;
	BED bedEntry;                                                                                                                        
	int lineNum = 0;

	// open the BED file for reading
	if (bedA->bedFile != "stdin") {

		ifstream bed(bedA->bedFile.c_str(), ios::in);
		if ( !bed ) {
			cerr << "Error: The requested bed file (" <<bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		
		BED a;
		while (getline(bed, bedLine)) {
			
			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			if (bedA->parseBedLine(a, bedFields, lineNum)) {

				vector<BED> hits;
				bedB->binKeeperFind(bedB->bedMap[a.chrom], a.start, a.end, hits);

				int numOverlaps = 0;
				for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {
					int s = max(a.start, h->start);
					int e = min(a.end, h->end);
					if (s < e) {
						if ( ((float)(e-s) / (float)(a.end - a.start)) > this->overlapFraction ) { 
							numOverlaps++;	
							if (!anyHit && !noHit && !writeCount) {			
								if (!writeB) {
									reportAIntersect(a,s,e);  cout << endl;
								}
								else {
									reportAIntersect(a,s,e); cout << "\t";
									reportB(*h); cout << endl;
								}
							}
						}
					}
				}
				if (anyHit && (numOverlaps >= 1)) {
					reportA(a); cout << endl;
				}
				else if (writeCount) {
					reportA(a); cout << "\t" << numOverlaps << endl;
				}
				else if (noHit && (numOverlaps == 0)) {
					reportA(a); cout << endl;
				}	
			}
		}
	}
	// "A" is being passed via STDIN.
	else {
		
		BED a;
		while (getline(cin, bedLine)) {

			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			if (bedA->parseBedLine(a, bedFields, lineNum)) {

				vector<BED> hits;
				bedB->binKeeperFind(bedB->bedMap[a.chrom], a.start, a.end, hits);

				int numOverlaps = 0;
				for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {

					int s = max(a.start, h->start);
					int e = min(a.end, h->end);
					if (s < e) {
						if ( ((float)(e-s) / (float)(a.end - a.start)) > this->overlapFraction ) { 
							numOverlaps++;	
							if (!anyHit && !noHit && !writeCount) {			
								if (!writeB) {
									reportAIntersect(a,s,e);  cout << endl;
								}
								else {
									reportAIntersect(a,s,e); cout << "\t";
									reportB(*h); cout << endl;
								}
							}
						} 
					}
					if (anyHit && (numOverlaps >= 1)) {
						reportA(a); cout << endl;
					}
					else if (writeCount) {
						reportA(a); cout << "\t" << numOverlaps << endl;
					}
					else if (noHit && (numOverlaps == 0)) {
						reportA(a); cout << endl;
					}
				}
			}
		}	
	}
}
// END Intersect BED3


