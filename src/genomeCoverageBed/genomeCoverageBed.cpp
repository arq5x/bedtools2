// 
//  genomeCoverageBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//
#include "lineFileUtilities.h"
#include "genomeCoverageBed.h"

/*
	Constructor
*/
BedCoverage::BedCoverage(string &bedFile, string &genomeFile, bool &eachBase, bool &startSites, int &max) {

	this->bedFile = bedFile;
	this->genomeFile = genomeFile;
	this->eachBase = eachBase;
	this->startSites = startSites;
	this->max = max;

	this->bed = new BedFile(bedFile);
}




/*
	Destructor
*/
BedCoverage::~BedCoverage(void) {

}




/*
	_Method: CoverageBeds
	_Purpose:
	_Example:
		           Chrom: ======================================
		 		   Reads:    -----	 -----
									-----		----
				   Depth: 00011112111111111001111000000000000000
	_Gotchas:
*/
void BedCoverage::CoverageBeds() {

	// Variables
	string chrom;
	unsigned int size;
	unsigned int genomeSize = 0;
	map<string, int, less<string> > chromSizes; 
	chromHistMap chromDepthHist;

	// open the GENOME file for reading.
	// if successful, load each chrom and it's size into
	// the "chromSize" map.  also compute the total size of the genom
	// and store in "genomeSize".
	ifstream genome(this->genomeFile.c_str(), ios::in);
	if ( !genome ) {
		cerr << "Error: The requested genome file (" <<this->genomeFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	else {
		while (genome >> chrom >> size) {
			chromSizes[chrom] = size;
			genomeSize += size;
		}
	}

	// open the BED file for reading.
	// if successful, store the start and end positions
	// for each entry (per chromosome)
	ifstream beds(this->bedFile.c_str(), ios::in);
	if ( !beds ) {
		cerr << "Error: The requested bed file (" <<this->bedFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	
	BED bedEntry;     // used to store the current BED line from the BED file.
	int lineNum = 0;
	string bedLine;	  // used to store the current (unparsed) line from the BED file.
	while (getline(beds, bedLine)) {
		
		vector<string> bedFields;
		Tokenize(bedLine,bedFields);

		lineNum++;
	
		if (this->bed->parseBedLine(bedEntry, bedFields, lineNum)) {
			this->chromCov[bedEntry.chrom][bedEntry.start+1].starts++;
			this->chromCov[bedEntry.chrom][bedEntry.end].ends++;
		}
	}

	// loop through each chromosome in the genome file "chromSizes".
	for (map<string, int, less<string> >::iterator m = chromSizes.begin(); m != chromSizes.end(); ++m) {
		
		// store the chromosome we are currently considering
		string chrom = m->first;

		if (this->chromCov.find(chrom) != this->chromCov.end()) {

			// does the user want to report the depth of each base?
			if (this->eachBase) {
				int depth = 0;  // initialize the depth
				
				for (int pos = 1; pos < chromSizes[chrom]; pos++) {

					// is there an entry in chromCov for this position?
					// if so, adjust the depth at the position by the number
					// of feature starts and ends
					if (this->chromCov[chrom].find(pos) != this->chromCov[chrom].end()) {
						depth += this->chromCov[chrom][pos].starts;
						
						// report the depth for this position.
						cout << chrom << "\t" << pos << "\t" << depth << endl;
						
						depth = depth - this->chromCov[chrom][pos].ends;
					}
					else {
						cout << chrom << "\t" << pos << "\t" << depth << endl;
					}
				}
			}
			// the user wants to create a histogram of coverage.
			else {
				int depth = 0;  // initialize the depth
				
				for (int pos = 1; pos < chromSizes[chrom]; pos++) {

					// is there an entry in chromCov for this position?
					// if so, adjust the depth at the position by the number
					// of feature starts and ends
					if (this->chromCov[chrom].find(pos) != this->chromCov[chrom].end()) {

						depth += this->chromCov[chrom][pos].starts;
					
						// add the depth at this position to the depth histogram
						// for this chromosome.  if the depth is greater than the
						// maximum bin requested, then readjust the depth to be the max
						if (depth >= this->max) {
							chromDepthHist[chrom][this->max]++;
						}
						else {
							chromDepthHist[chrom][depth]++;
						}
					}
					// otherwise, just track the current (unadjusted) depth 
					// OR the max if depth >= max.
					else {
						if (depth >= this->max) {
							chromDepthHist[chrom][this->max]++;
						}
						else {
							chromDepthHist[chrom][depth]++;
						}
					}
					depth = depth - this->chromCov[chrom][pos].ends;
				}
				
				// report the histogram for each chromosome
				for (histMap::iterator depthIt = chromDepthHist[chrom].begin(); depthIt != chromDepthHist[chrom].end(); ++depthIt) {
					int depth = depthIt->first;
					unsigned int numBasesAtDepth = depthIt->second;  
					
					cout << chrom << "\t" << depth << "\t" << numBasesAtDepth << "\t" 
						<< chromSizes[chrom] << "\t" << (float) ((float)numBasesAtDepth / (float)chromSizes[chrom]) << endl;
				}
			}
		}
		// There are not BED entries for this chromosome.
		// We therefore set the depth for each base therein to 0
		else {
			chromDepthHist[chrom][0] += chromSizes[chrom];
			
			// report the histogram for each chromosome
			for (histMap::iterator depthIt = chromDepthHist[chrom].begin(); depthIt != chromDepthHist[chrom].end(); ++depthIt) {
				int depth = depthIt->first;
				unsigned int numBasesAtDepth = depthIt->second;
				
				cout << chrom << "\t" << depth << "\t" << numBasesAtDepth << "\t" 
					<< chromSizes[chrom] << "\t" << (float) ((float)numBasesAtDepth / (float)chromSizes[chrom]) << endl;
			}
		}
	}
	
	// report the histogram for the genome as a whole: assuming the user
	// has not requested "per base coverage".
	if (!this->eachBase) {
		
		histMap genomeHist;  // depth histogram for the entire genome
		
		// loop through each chromosome and add the depth and number of bases at each depth
		// to the aggregate histogram for the entire genome
		for (chromHistMap::iterator chromIt = chromDepthHist.begin(); chromIt != chromDepthHist.end(); ++chromIt) {
			string chrom = chromIt->first;
			for (histMap::iterator depthIt = chromDepthHist[chrom].begin(); depthIt != chromDepthHist[chrom].end(); ++depthIt) {
				int depth = depthIt->first;
				unsigned int numBasesAtDepth = depthIt->second;
				
				genomeHist[depth] += numBasesAtDepth;
			}
		}
		
		// loop through the depths for the entire genome
		// and report the number and fraction of bases in
		// the entire genome that are at said depth.
		for (histMap::iterator genomeDepthIt = genomeHist.begin(); genomeDepthIt != genomeHist.end(); ++genomeDepthIt) {
			int depth = genomeDepthIt->first;
			unsigned int numBasesAtDepth = genomeDepthIt->second;
			
			cout << "genome" << "\t" << depth << "\t" << numBasesAtDepth << "\t" 
				<< genomeSize << "\t" << (float) ((float)numBasesAtDepth / (float)genomeSize) << endl;
		}	
	}
}

