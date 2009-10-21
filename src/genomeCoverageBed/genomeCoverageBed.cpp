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
		string chrom;
		unsigned int size;
		while (genome >> chrom >> size) {
			chromSizes[chrom] = size;
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
	string prevChrom, currChrom;
	vector<DEPTH> chromCov;
	int chromSize = 0;
	int start, end;
	
	while (getline(beds, bedLine)) {
		
		vector<string> bedFields;
		Tokenize(bedLine,bedFields);
		lineNum++;

		if (this->bed->parseBedLine(bedEntry, bedFields, lineNum)) {
						
			currChrom = bedEntry.chrom;
			start = bedEntry.start;
			end = bedEntry.end -1;
			
			if (currChrom != prevChrom)  {
				// If we've moved beyond the first encountered chromosomes,
				// process the results of the previous chromosome.
				if (prevChrom.length() > 0) {
					ReportChromCoverage(chromCov, chromSizes[prevChrom], prevChrom, chromDepthHist);
				}
				
				// empty the previous chromosome and reserve new
				// space for the current one
				std::vector< DEPTH >().swap(chromCov);
				chromSize = chromSizes[currChrom];
				chromCov.resize(chromSize);

				// process the first line for this chromosome.
				// make sure the coordinates fit within the chrom
				if (start < chromSize) {
					chromCov[start].starts++;
				}
				if (end < chromSize) {
					chromCov[end].ends++;
				}
				else {
					chromCov[chromSize-1].ends++;
				}
			}
			else {
				// process the other lines for this chromosome.
				// make sure the coordinates fit within the chrom
				if (start < chromSize) {
					chromCov[start].starts++;
				}
				if (end < chromSize) {
					chromCov[end].ends++;
				}
				else {
					chromCov[chromSize-1].ends++;
				}			
			}
			prevChrom = currChrom;
		}
	}
	// process the results of the last chromosome.
	ReportChromCoverage(chromCov, chromSizes[currChrom], currChrom, chromDepthHist);
	
	if (! this->eachBase) {
		ReportGenomeCoverage(chromSizes, chromDepthHist);
	}	
}



void BedCoverage::ReportChromCoverage(vector<DEPTH> &chromCov, int &chromSize, string &chrom, chromHistMap &chromDepthHist) {
	
	if (this->eachBase) {
		int depth = 0;  // initialize the depth
		for (int pos = 0; pos < chromSize; pos++) {
			depth += chromCov[pos].starts;
			
			// report the depth for this position.
			cout << chrom << "\t" << pos+1 << "\t" << depth << endl;
			
			depth = depth - chromCov[pos].ends;
		}
	}
	else {
		
		int depth = 0;  // initialize the depth
		
		for (int pos = 0; pos < chromSize; pos++) {
			
			depth += chromCov[pos].starts;
			
			// add the depth at this position to the depth histogram
			// for this chromosome.  if the depth is greater than the
			// maximum bin requested, then readjust the depth to be the max
			if (depth >= this->max) {
				chromDepthHist[chrom][this->max]++;
			}
			else {
				chromDepthHist[chrom][depth]++;
			}
			depth = depth - chromCov[pos].ends;
		}
		// report the histogram for each chromosome
		for (histMap::iterator depthIt = chromDepthHist[chrom].begin(); depthIt != chromDepthHist[chrom].end(); ++depthIt) {
			int depth = depthIt->first;
			unsigned int numBasesAtDepth = depthIt->second;  
			
			cout << chrom << "\t" << depth << "\t" << numBasesAtDepth << "\t" 
				<< chromSize << "\t" << (float) ((float)numBasesAtDepth / (float)chromSize) << endl;
		}
	}
}



void BedCoverage::ReportGenomeCoverage(map<string, int> &chromSizes, chromHistMap &chromDepthHist) {
	
	unsigned int genomeSize = 0;
	for (map<string, int, less<string> >::iterator m = chromSizes.begin(); m != chromSizes.end(); ++m) {
		
		string chrom = m->first;
		genomeSize += chromSizes[chrom];
		// if there were no reads for a give chromosome, then
		// add the length of the chrom to the 0 bin.
		if ( chromDepthHist.find(chrom) == chromDepthHist.end() ) {
			chromDepthHist[chrom][0] += chromSizes[chrom];
		}
	}

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
