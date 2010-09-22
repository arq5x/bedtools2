/*****************************************************************************
  complementBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "complementBed.h"

BedComplement::BedComplement(string &bedFile, string &genomeFile) {

	_bedFile = bedFile;
	_genomeFile = genomeFile;
	
	_bed    = new BedFile(bedFile);
	_genome = new GenomeFile(genomeFile);
		
}


BedComplement::~BedComplement(void) {
}


//
// Merge overlapping BED entries into a single entry 
//
void BedComplement::ComplementBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	_bed->loadBedFileIntoMapNoBin();
	
	vector<short> chromMasks;
	string currChrom;
	
	// loop through each chromosome and merge their BED entries
	masterBedMapNoBin::const_iterator m    = _bed->bedMapNoBin.begin();
	masterBedMapNoBin::const_iterator mEnd = _bed->bedMapNoBin.end();
    for (; m != mEnd; ++m) {
		currChrom = m->first;
		CHRPOS currChromSize = _genome->getChromSize(currChrom);
		
		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 
		
		// create a flag for every base on the chrom.
		vector<short> chromMasks(currChromSize, 0);
		
		vector<BED>::const_iterator bIt  = bedList.begin();
		vector<BED>::const_iterator bEnd = bedList.end();
		for ( ; bIt != bEnd; ++bIt) {
			
			// sanity check the end of the bed entry
			if (bIt->end > currChromSize) {
				cout << "End of BED entry exceeds chromosome length. Please correct." << endl;
				_bed->reportBedNewLine(*bIt);
				exit(1);
			}
			
			// mask all of the positions spanned by this BED entry.
			for (CHRPOS b = bIt->start; b < bIt->end; b++)
				chromMasks[b] = 1;
		}
		
		CHRPOS i = 0;
		CHRPOS start;
		while (i < chromMasks.size()) {
			if (chromMasks[i] == 0) {
				start = i;
				while ((chromMasks[i] == 0) && (i < chromMasks.size()))
					i++;
				
				if (start > 0) 
				    cout << currChrom << "\t" << start << "\t" << i << endl;
				else 
				    cout << currChrom << "\t" << 0 << "\t" << i << endl;
			}
			i++;
		}
	}
}

