/*****************************************************************************
  genomeCoverage.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "genomeCoverageBed.h"

/*
	Constructor
*/
BedCoverage::BedCoverage(string &bedFile, string &genomeFile, bool &eachBase, bool &startSites, bool &bedGraph, int &max) {

	this->bedFile    = bedFile;
	this->genomeFile = genomeFile;
	this->eachBase   = eachBase;
	this->startSites = startSites;
	this->bedGraph   = bedGraph;
	this->max        = max;

	this->bed        = new BedFile(bedFile);
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
void BedCoverage::CoverageBeds(istream &bedInput) {

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

	string prevChrom, currChrom;
	vector<DEPTH> chromCov;
	int chromSize = 0;
	int start, end;
	
	string bedLine;                                                                                                                    
	int lineNum = 0;					// current input line number
	vector<string> bedFields;			// vector for a BED entry
	bedFields.reserve(12);
	
	while (getline(bedInput, bedLine)) {
		
		Tokenize(bedLine,bedFields);
		lineNum++;
		BED bedEntry;     // used to store the current BED line from the BED file.
		
		if (bed->parseLine(bedEntry, bedFields, lineNum)) {
						
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
		bedFields.clear();
	}
	// process the results of the last chromosome.
	ReportChromCoverage(chromCov, chromSizes[currChrom], currChrom, chromDepthHist);
	
	if (this->eachBase == false && this->bedGraph == false) {
		ReportGenomeCoverage(chromSizes, chromDepthHist);
	}
}



void BedCoverage::ReportChromCoverage(const vector<DEPTH> &chromCov, int &chromSize, string &chrom, chromHistMap &chromDepthHist) {
	
	if (this->eachBase) {
		int depth = 0;  // initialize the depth
		for (int pos = 0; pos < chromSize; pos++) {
			depth += chromCov[pos].starts;
			
			// report the depth for this position.
			cout << chrom << "\t" << pos+1 << "\t" << depth << endl;
			
			depth = depth - chromCov[pos].ends;
		}
	}
	else if (this->bedGraph) {
		ReportChromCoverageBedGraph(chromCov, chromSize, chrom);
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


void BedCoverage::ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, int &chromSize, string &chrom) {

	int depth     = 0;     // initialize the depth
	int lastStart = -1 ;
	int lastDepth = -1 ;

	for (int pos = 0; pos < chromSize; pos++) {
		depth += chromCov[pos].starts;

		if (depth == 0 && lastDepth != -1) {
			// We've found a new block of zero coverage, so report
			// the previous block of non-zero coverage.
			cout << chrom << "\t" << lastStart << "\t" << pos << "\t" << lastDepth << endl;
			lastDepth = -1;
			lastStart = -1;
		}
		else if (depth > 0 && depth != lastDepth) {
			// Coverage depth has changed, print the last interval coverage (if any)
			if (lastDepth != -1) { 
				cout << chrom << "\t" << lastStart << "\t" << pos << "\t" << lastDepth << endl;
			}
			//Set current position as the new interval start + depth
			lastDepth = depth;
			lastStart = pos;
		}
		// Default: the depth has not changed, so we will not print anything.
		// Proceed until the depth changes.
		
		// Update depth
		depth = depth - chromCov[pos].ends;
	}
	
	
	//Print information about the last position
	if (lastDepth != -1) {
		cout << chrom << "\t" << lastStart << "\t" << chromSize << "\t" << lastDepth << endl;
	}
}


void BedCoverage::DetermineBedInput() {
	if (bed->bedFile != "stdin") {   // process a file
		ifstream beds(bed->bedFile.c_str(), ios::in);
		if ( !beds ) {
			cerr << "Error: The requested bed file (" << bed->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		CoverageBeds(beds);
	}
	else {   						// process stdin
		CoverageBeds(cin);		
	}
}
