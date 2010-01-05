/*****************************************************************************
  complementBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "complementBed.h"

//==========================
// Constructor
//
BedComplement::BedComplement(string &bedFile, string &genomeFile) {

	this->bedFile = bedFile;
	this->genomeFile = genomeFile;
	this->bed = new BedFile(bedFile);
}


//
// Destructor
//
BedComplement::~BedComplement(void) {
}


//
// Merge overlapping BED entries into a single entry 
//
void BedComplement::ComplementBed() {

	// open the GENOME file for reading
	ifstream genome(this->genomeFile.c_str(), ios::in);
	if ( !genome ) {
		cerr << "Error: The requested genome file (" <<this->genomeFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}

	string chrom;
	unsigned int size;

	map<string, int, less<string> > chromSizes; 

	while (genome >> chrom >> size) {
		chromSizes[chrom] = size;
	}

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();
	
	vector<short> chromMasks;
	string currChrom;
	
	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		currChrom = m->first;
		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 
		
		// create a flag for every base on the chrom.
		vector<short> chromMasks(chromSizes[currChrom], 0);
		
		vector<BED>::const_iterator bIt = bedList.begin();
		vector<BED>::const_iterator bEnd = bedList.end();
		for ( ; bIt != bEnd; ++bIt) {
			
			// sanity check the end of the bed entry
			if (bIt->end > chromSizes[currChrom]) {
				cout << "End of BED entry exceeds chromosome length. Please correct." << endl;
				bed->reportBedNewLine(*bIt);
				exit(1);
			}
			
			// mask all of the positions spanned by this BED entry.
			for (int b = bIt->start; b < bIt->end; b++) {
				chromMasks[b] = 1;
			}
		}
		
		unsigned int i = 0;
		unsigned int start;
		while (i < chromMasks.size()) {
			if (chromMasks[i] == 0) {
				start = i;
				while ((chromMasks[i] == 0) && (i < chromMasks.size())) {
					i++;
				}
				
				if (start > 0) {
					cout << currChrom << "\t" << start << "\t" << i << endl;
				}
				else {
					cout << currChrom << "\t" << 0 << "\t" << i << endl;
				}
			}
			i++;
		}
	}
}

