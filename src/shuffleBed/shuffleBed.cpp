/*****************************************************************************
  shuffleBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "shuffleBed.h"


BedShuffle::BedShuffle(string &bedFile, string &genomeFile, string &excludeFile, bool &haveSeed, bool &haveExclude, bool &sameChrom, int &seed) {

	this->bedFile     = bedFile;
	this->genomeFile  = genomeFile;
	this->excludeFile = excludeFile;
	this->sameChrom   = sameChrom;
	this->haveExclude = haveExclude;
	this->haveSeed    = haveSeed;

	
	// use the supplied seed for the random
	// number generation if given.  else,
	// roll our own.
	if (this->haveSeed) {
		this->seed = seed;
		srand(seed);
	}
	else {
		srand((unsigned)time(0)); 
	}
	
	this->bed         = new BedFile(bedFile);
	this->genome      = new GenomeFile(genomeFile);
	this->chroms      = genome->getChromList();
	this->numChroms   = genome->getNumberOfChroms();
	
	if (this->haveExclude) {
		this->exclude = new BedFile(excludeFile);
		this->exclude->loadBedFileIntoMap();	
	}	
}


BedShuffle::~BedShuffle(void) {

}


void BedShuffle::Shuffle(istream &bedInput) {

	int lineNum = 0;
	string bedLine;	  // used to store the current (unparsed) line from the BED file.
	vector<string> bedFields;
	bedFields.reserve(12);

	while (getline(bedInput, bedLine)) {
		
		Tokenize(bedLine,bedFields);
		lineNum++;
		BED bedEntry;     // used to store the current BED line from the BED file.
				
		if (bed->parseLine(bedEntry, bedFields, lineNum)) {
			// choose a new locus for this feat
			ChooseLocus(bedEntry);			
			bed->reportBedNewLine(bedEntry);
		}
		bedFields.clear();	
	}
}



void BedShuffle::ShuffleWithExclusions(istream &bedInput) {

	int lineNum = 0;
	string bedLine;	  // used to store the current (unparsed) line from the BED file.
	vector<string> bedFields;
	bedFields.reserve(12);		
	vector<BED> hits;
	hits.reserve(100);
		
	while (getline(bedInput, bedLine)) {
		
		Tokenize(bedLine,bedFields);
		lineNum++;
		BED bedEntry;     // used to store the current BED line from the BED file.

		if (bed->parseLine(bedEntry, bedFields, lineNum)) {
						
			// choose a random locus
			ChooseLocus(bedEntry);	
			
			// test to see if the chosen locus overlaps 
			// with an exclude region
			exclude->FindOverlapsPerBin(bedEntry.chrom, bedEntry.start, bedEntry.end, bedEntry.strand, hits, false);
					
			bool haveOverlap = false;
			vector<BED>::const_iterator hitsItr = hits.begin();
			vector<BED>::const_iterator hitsEnd = hits.end();
			for (; hitsItr != hitsEnd; ++hitsItr) {

				int s = max(bedEntry.start, hitsItr->start);
				int e = min(bedEntry.end, hitsItr->end);

				if ( (e - s) > 0) {
					haveOverlap = true;
					break;   /* stop looking.  one overlap is enough*/
				}
			}
			
			/* 
			   keep looking as long as the chosen
			   locus happens to overlap with regions
			   that the user wishes to exclude.
			*/
			int tries = 0;
			while ((haveOverlap == true) && (tries <= MAX_TRIES)) {

				// choose a new locus
				ChooseLocus(bedEntry);

				vector<BED> hits;
				exclude->FindOverlapsPerBin(bedEntry.chrom, bedEntry.start, bedEntry.end, 
											bedEntry.strand, hits, false);

				haveOverlap = false;
				vector<BED>::const_iterator hitsItr = hits.begin();
				vector<BED>::const_iterator hitsEnd = hits.end();
				for (; hitsItr != hitsEnd; ++hitsItr) {

					int s = max(bedEntry.start, hitsItr->start);
					int e = min(bedEntry.end, hitsItr->end);

					if ( (e - s) > 0) {
						haveOverlap = true;
						break;  /* stop looking.  one overlap is enough*/
					}
				}
				tries++;
			}
			
			if (tries > MAX_TRIES) {
				cerr << "Error, line " << lineNum << ": tried " << MAX_TRIES << " potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on." << endl;
			}
			else {
				bed->reportBedNewLine(bedEntry);
			}
		}
		bedFields.clear();
	}
}



void BedShuffle::ChooseLocus(BED &bedEntry) {
	
	string chrom = bedEntry.chrom;
	int start    = bedEntry.start;
	int end      = bedEntry.end;
	int length   = end - start;
	
	string randomChrom;
	int randomStart;
	int chromSize;
	
	if (!this->sameChrom) {
		randomChrom    = chroms[rand() % this->numChroms];
		chromSize      = genome->getChromSize(randomChrom);
		randomStart    = rand() % chromSize;
		bedEntry.chrom = randomChrom;
		bedEntry.start = randomStart;
		bedEntry.end   = randomStart + length;
	}
	else {
		chromSize      = genome->getChromSize(chrom);
		randomStart    = rand() % chromSize;
		bedEntry.start = randomStart;
		bedEntry.end   = randomStart + length;
	}
	
	// ensure that the chosen location doesn't go past
	// the length of the chromosome. if so, keep looking
	// for a new spot.
	while (bedEntry.end > chromSize) {
		if (!this->sameChrom) {
			randomChrom    = chroms[rand() % this->numChroms];
			chromSize      = genome->getChromSize(randomChrom);
			randomStart    = rand() % chromSize;
			bedEntry.chrom = randomChrom;
			bedEntry.start = randomStart;
			bedEntry.end   = randomStart + length;
		}
		else {
			chromSize      = genome->getChromSize(chrom);
			randomStart    = rand() % chromSize;
			bedEntry.start = randomStart;
			bedEntry.end   = randomStart + length;
		}
	}
}


void BedShuffle::DetermineBedInput() {
	if (bed->bedFile != "stdin") {   // process a file
		ifstream beds(bed->bedFile.c_str(), ios::in);
		if ( !beds ) {
			cerr << "Error: The requested bed file (" << bed->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		if (this->haveExclude) { 
			ShuffleWithExclusions(beds); 
		}
		else {
			Shuffle(beds); 
		}
	}
	else {
									// process stdin
		if (this->haveExclude) { 
			ShuffleWithExclusions(cin); 
		}
		else {
			Shuffle(cin); 
		}		
	}
}

