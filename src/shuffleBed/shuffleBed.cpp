/*****************************************************************************
  shuffleBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "shuffleBed.h"


BedShuffle::BedShuffle(string &bedFile, string &genomeFile, string &excludeFile, bool &haveSeed, bool &haveExclude, bool &sameChrom, int &seed) {

	_bedFile     = bedFile;
	_genomeFile  = genomeFile;
	_excludeFile = excludeFile;
	_sameChrom   = sameChrom;
	_haveExclude = haveExclude;
	_haveSeed    = haveSeed;

	
	// use the supplied seed for the random
	// number generation if given.  else,
	// roll our own.
	if (_haveSeed) {
		_seed = seed;
		srand(seed);
	}
	else {
		srand((unsigned)time(0)); 
	}
	
	_bed         = new BedFile(bedFile);
	_genome      = new GenomeFile(genomeFile);
	_chroms      = _genome->getChromList();
	_numChroms   = _genome->getNumberOfChroms();
	
	if (_haveExclude) {
		_exclude = new BedFile(excludeFile);
		_exclude->loadBedFileIntoMap();	
	}
	
	if (_bed->bedFile != "stdin") {   // process a file
		if (_haveExclude)
			ShuffleWithExclusions(); 
		else
			Shuffle(); 
	}
	else {				// process stdin
		if (_haveExclude)
			ShuffleWithExclusions(); 
		else
			Shuffle(); 
	}	
}


BedShuffle::~BedShuffle(void) {

}


void BedShuffle::Shuffle() {

	int lineNum = 0;
	BED bedEntry, nullBed;     // used to store the current BED line from the BED file.
	BedLineStatus bedStatus;
	
	_bed->Open();
	while ((bedStatus = _bed->GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
		if (bedStatus == BED_VALID) {
			ChooseLocus(bedEntry);			
			_bed->reportBedNewLine(bedEntry);
			bedEntry = nullBed;
		}
	}
	_bed->Close();
}



void BedShuffle::ShuffleWithExclusions() {

	int lineNum = 0;
	BED bedEntry, nullBed;     // used to store the current BED line from the BED file.
	BedLineStatus bedStatus;
	vector<BED> hits;
	hits.reserve(100);
		
	_bed->Open();	
	while ((bedStatus = _bed->GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
		if (bedStatus == BED_VALID) {				
    		// choose a random locus
    		ChooseLocus(bedEntry);	
		
    		// test to see if the chosen locus overlaps 
    		// with an exclude region
    		_exclude->FindOverlapsPerBin(bedEntry.chrom, bedEntry.start, bedEntry.end, bedEntry.strand, hits, false);
				
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
    			_exclude->FindOverlapsPerBin(bedEntry.chrom, bedEntry.start, bedEntry.end, 
    										bedEntry.strand, hits, false);

    			haveOverlap = false;
    			vector<BED>::const_iterator hitsItr = hits.begin();
    			vector<BED>::const_iterator hitsEnd = hits.end();
    			for (; hitsItr != hitsEnd; ++hitsItr) {

    				int s = max(bedEntry.start, hitsItr->start);
    				int e = min(bedEntry.end, hitsItr->end);

    				if ( (e - s) > 0) {
    					haveOverlap = true;
    					break;  // stop looking.  one overlap is enough
    				}
    			}
    			tries++;
    		}
		
    		if (tries > MAX_TRIES) {
    			cerr << "Error, line " << lineNum << ": tried " << MAX_TRIES << " potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on." << endl;
    		}
    		else {
    			_bed->reportBedNewLine(bedEntry);
    		}
    	}
    	bedEntry = nullBed;
	}
	_bed->Close();
}



void BedShuffle::ChooseLocus(BED &bedEntry) {
	
	string chrom = bedEntry.chrom;
	CHRPOS start    = bedEntry.start;
	CHRPOS end      = bedEntry.end;
	CHRPOS length   = end - start;
	
	string randomChrom;
	CHRPOS randomStart;
	CHRPOS chromSize;
	
	if (_sameChrom == false) {
		randomChrom    = _chroms[rand() % _numChroms];
		chromSize      = _genome->getChromSize(randomChrom);
		randomStart    = rand() % chromSize;
		bedEntry.chrom = randomChrom;
		bedEntry.start = randomStart;
		bedEntry.end   = randomStart + length;
	}
	else {
		chromSize      = _genome->getChromSize(chrom);
		randomStart    = rand() % chromSize;
		bedEntry.start = randomStart;
		bedEntry.end   = randomStart + length;
	}
	
	// ensure that the chosen location doesn't go past
	// the length of the chromosome. if so, keep looking
	// for a new spot.
	while (bedEntry.end > chromSize) {
		if (_sameChrom == false) {
			randomChrom    = _chroms[rand() % _numChroms];
			chromSize      = _genome->getChromSize(randomChrom);
			randomStart    = rand() % chromSize;
			bedEntry.chrom = randomChrom;
			bedEntry.start = randomStart;
			bedEntry.end   = randomStart + length;
		}
		else {
			chromSize      = _genome->getChromSize(chrom);
			randomStart    = rand() % chromSize;
			bedEntry.start = randomStart;
			bedEntry.end   = randomStart + length;
		}
	}
}

