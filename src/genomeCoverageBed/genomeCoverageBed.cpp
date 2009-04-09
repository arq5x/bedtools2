#include "genomeCoverageBed.h"

BedCoverage::BedCoverage(string &bedFile, string &genomeFile, bool &eachBase, bool &startSites, int &max) {

	this->bedFile = bedFile;
	this->genomeFile = genomeFile;
	this->eachBase = eachBase;
	this->startSites = startSites;
	this->max = max;

	this->bed = new BedFile(bedFile);
}


BedCoverage::~BedCoverage(void) {
}

void BedCoverage::CoverageBeds() {

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

	bed->loadBedFileIntoMapNoBin();
	for (map<string, int, less<string> >::iterator m = chromSizes.begin(); m != chromSizes.end(); ++m) {

		if (bed->bedMapNoBin.find(m->first) != bed->bedMapNoBin.end()) {

			// declare a vector with an element for each bp on the current chromosome
			vector<COV> chromCov(chromSizes[m->first]);

			// get the BED entries for this chromosome
			vector<BED> bedList = bed->bedMapNoBin[m->first];

			// record the start and end sites of each bed entry
			for (int i = 0; i < bedList.size(); ++i) {			
				chromCov[bedList[i].start - 1].starts++;
				chromCov[bedList[i].end - 1].ends++;
			}

			// does the user want to report the depth of each base?
			if (this->eachBase) {
				int d = 0;
				for (int pos = 0; pos < chromCov.size(); pos++) {

					d += (chromCov[pos].starts - chromCov[pos].ends);
					chromCov[pos].depth = d + chromCov[pos].ends;

					cout << m->first << "\t" << pos+1 << "\t" << chromCov[pos].depth << endl; 
				}	
			}
			else {
				int d = 0;
				map<int, int, less<int> > depthHist;
				for (int pos = 0; pos < chromCov.size(); pos++) {

					d += (chromCov[pos].starts - chromCov[pos].ends);
					chromCov[pos].depth = d + chromCov[pos].ends;

					if ((d + chromCov[pos].ends) >= this->max)
						depthHist[this->max]++;
					else {
						depthHist[d + chromCov[pos].ends]++;
					}
				}

				// report the histogram
				for (map<int, int, less<string> >::iterator dI = depthHist.begin(); dI != depthHist.end(); ++dI) {
					cout << m->first << "\t" << dI->first << "\t" << dI->second << "\t" 
						<< chromSizes[m->first] << "\t" << (float) ((float)dI->second / (float)chromSizes[m->first]) << endl;
				}
			}
		}
	}
}

