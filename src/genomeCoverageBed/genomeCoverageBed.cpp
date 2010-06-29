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


BedGenomeCoverage::BedGenomeCoverage(string bedFile, string genomeFile, bool eachBase, 
	                                 bool startSites, bool bedGraph, bool bedGraphAll, int max, bool bamInput) {

	_bedFile       = bedFile;
	_genomeFile    = genomeFile;
	_eachBase      = eachBase;
	_startSites    = startSites;
	_bedGraph      = bedGraph;
	_bedGraphAll   = bedGraphAll;
	_max           = max;
	_bamInput      = bamInput;
	
	_bed        = new BedFile(bedFile);
	_genome     = new GenomeFile(genomeFile);
	
	if (_bed->bedFile != "stdin") {   // process a file
		if (_bamInput == false)
			CoverageBed();
		else 
			CoverageBam(_bed->bedFile);
	}
	else {   // process stdin
		if (_bamInput == false) 
			CoverageBed();
		else 
			CoverageBam("stdin");	
	}
}


BedGenomeCoverage::~BedGenomeCoverage(void) {
	delete _bed;
	delete _genome;
}


void BedGenomeCoverage::CoverageBed() {

	chromHistMap chromDepthHist;

	string prevChrom, currChrom;
	vector<DEPTH> chromCov;

	int prevChromSize = 0;
	int currChromSize = 0;
	int start, end;
	
	BED a, nullBed;
	int lineNum = 0;					// current input line number
	BedLineStatus bedStatus;
	
	_bed->Open();
	bedStatus = _bed->GetNextBed(a, lineNum);
	while (bedStatus != BED_INVALID) {
		// ignore headers and blank lines
		if (bedStatus == BED_VALID) {	
			currChrom = a.chrom;
			start     = a.start;
			end       = a.end - 1;
				
			if (currChrom != prevChrom)  {
				// If we've moved beyond the first encountered chromosomes,
				// process the results of the previous chromosome.
				if (prevChrom.length() > 0) {
					ReportChromCoverage(chromCov, prevChromSize, prevChrom, chromDepthHist);
				}
			
				// empty the previous chromosome and reserve new
				std::vector<DEPTH>().swap(chromCov);
			
				// get the current chrom size and allocate space 
				currChromSize = _genome->getChromSize(currChrom);
				if (currChromSize >= 0)
					chromCov.resize(currChromSize);
				else {
					cerr << "Chromosome " << currChrom << " found in your BED file but not in your genome file.  Exiting." << endl;
					exit(1);
				}
					
				// process the first line for this chromosome.
				// make sure the coordinates fit within the chrom
				if (start < currChromSize) {
					chromCov[start].starts++;
				}
				if (end < currChromSize) {
					chromCov[end].ends++;
				}
				else {
					chromCov[currChromSize-1].ends++;
				}
			}
			else {
				// process the other lines for this chromosome.
				// make sure the coordinates fit within the chrom
				if (start < currChromSize) {
					chromCov[start].starts++;
				}
				if (end < currChromSize) {
					chromCov[end].ends++;
				}
				else {
					chromCov[currChromSize-1].ends++;
				}			
			}
			prevChrom     = currChrom;
			prevChromSize = currChromSize;
			a = nullBed;
		}
		bedStatus = _bed->GetNextBed(a, lineNum);
	}
	_bed->Close();

	// process the results of the last chromosome.
	ReportChromCoverage(chromCov, currChromSize, currChrom, chromDepthHist);	
	if (_eachBase == false && _bedGraph == false && _bedGraphAll == false) {
		ReportGenomeCoverage(chromDepthHist);
	}
}


void BedGenomeCoverage::CoverageBam(string bamFile) {

	chromHistMap chromDepthHist;

	string prevChrom, currChrom;
	vector<DEPTH> chromCov;
	
	int prevChromSize = 0;
	int currChromSize = 0;
	int start, end;
	
	// open the BAM file
	BamReader reader;
	reader.Open(bamFile);

	// get header & reference information
	string header = reader.GetHeaderText();
	RefVector refs = reader.GetReferenceData();

	// convert each aligned BAM entry to BED 
	// and compute coverage on B
	BamAlignment bam;	
	while (reader.GetNextAlignment(bam)) {
		
		if (bam.IsMapped()) {
			
			currChrom  = refs.at(bam.RefID).RefName;
			start      = bam.Position;
			end        = bam.GetEndPosition(false) - 1;
			
			if (currChrom != prevChrom)  {
				// If we've moved beyond the first encountered chromosomes,
				// process the results of the previous chromosome.
				if (prevChrom.length() > 0) {
					ReportChromCoverage(chromCov, prevChromSize, prevChrom, chromDepthHist);
				}
				
				// empty the previous chromosome and reserve new
				std::vector<DEPTH>().swap(chromCov);

				// get the current chrom size and allocate space 
				currChromSize = _genome->getChromSize(currChrom);				
				if (currChromSize >= 0)
					chromCov.resize(currChromSize);
				else {
					cerr << "Chromosome " << currChrom << " found in your BED file but not in your genome file.  Exiting." << endl;
					exit(1);
				}

				// process the first line for this chromosome.
				// make sure the coordinates fit within the chrom
				if (start < currChromSize) {
					chromCov[start].starts++;
				}
				if (end < currChromSize) {
					chromCov[end].ends++;
				}
				else {
					chromCov[currChromSize-1].ends++;
				}
			}
			else {
				// process the other lines for this chromosome.
				// make sure the coordinates fit within the chrom
				if (start < currChromSize) {
					chromCov[start].starts++;
				}
				if (end < currChromSize) {
					chromCov[end].ends++;
				}
				else {
					chromCov[currChromSize-1].ends++;
				}			
			}
			prevChrom     = currChrom;
			prevChromSize = currChromSize;
		}
	}
	// process the results of the last chromosome.
	ReportChromCoverage(chromCov, currChromSize, currChrom, chromDepthHist);
	
	if (_eachBase == false && _bedGraph == false && _bedGraphAll == false) {
		ReportGenomeCoverage(chromDepthHist);
	}
	
	// close the BAM
	reader.Close();
}


void BedGenomeCoverage::ReportChromCoverage(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom, chromHistMap &chromDepthHist) {
	
	if (_eachBase) {
		int depth = 0;  // initialize the depth
		for (int pos = 0; pos < chromSize; pos++) {
			
			depth += chromCov[pos].starts;
			// report the depth for this position.
			cout << chrom << "\t" << pos+1 << "\t" << depth << endl;
			depth = depth - chromCov[pos].ends;
		}
	}
	else if (_bedGraph == true || _bedGraphAll == true) {
		ReportChromCoverageBedGraph(chromCov, chromSize, chrom);
	}
	else {
		
		int depth = 0;  // initialize the depth
		
		for (int pos = 0; pos < chromSize; pos++) {
			
			depth += chromCov[pos].starts;
			
			// add the depth at this position to the depth histogram
			// for this chromosome.  if the depth is greater than the
			// maximum bin requested, then readjust the depth to be the max
			if (depth >= _max) {
				chromDepthHist[chrom][_max]++;
			}
			else {
				chromDepthHist[chrom][depth]++;
			}
			depth = depth - chromCov[pos].ends;
		}
		// report the histogram for each chromosome
		histMap::const_iterator depthIt  = chromDepthHist[chrom].begin();
		histMap::const_iterator depthEnd = chromDepthHist[chrom].end();
		for (; depthIt != depthEnd; ++depthIt) {
			int depth                    = depthIt->first;
			unsigned int numBasesAtDepth = depthIt->second;  
			cout << chrom << "\t" << depth << "\t" << numBasesAtDepth << "\t" 
				<< chromSize << "\t" << (float) ((float)numBasesAtDepth / (float)chromSize) << endl;
		}
	}
}



void BedGenomeCoverage::ReportGenomeCoverage(chromHistMap &chromDepthHist) {
	
	// get the list of chromosome names in the genome
	vector<string> chromList = _genome->getChromList();

	unsigned int genomeSize = 0;
	vector<string>::const_iterator chromItr = chromList.begin();
	vector<string>::const_iterator chromEnd = chromList.end();	
	for (; chromItr != chromEnd; ++chromItr) {	
		string chrom = *chromItr;
		genomeSize   += _genome->getChromSize(chrom);
		// if there were no reads for a give chromosome, then
		// add the length of the chrom to the 0 bin.
		if ( chromDepthHist.find(chrom) == chromDepthHist.end() ) {
			chromDepthHist[chrom][0] += _genome->getChromSize(chrom);
		}
	}

	histMap genomeHist;  // depth histogram for the entire genome
	
	// loop through each chromosome and add the depth and number of bases at each depth
	// to the aggregate histogram for the entire genome
	for (chromHistMap::iterator chromIt = chromDepthHist.begin(); chromIt != chromDepthHist.end(); ++chromIt) {
		string chrom = chromIt->first;
		for (histMap::iterator depthIt = chromDepthHist[chrom].begin(); depthIt != chromDepthHist[chrom].end(); ++depthIt) {
			int depth                    = depthIt->first;
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


void BedGenomeCoverage::ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom) {

	int depth     = 0;     // initialize the depth
	int lastStart = -1;
	int lastDepth = -1;

	for (int pos = 0; pos < chromSize; pos++) {
		depth += chromCov[pos].starts;

		if (depth != lastDepth) {
			// Coverage depth has changed, print the last interval coverage (if any)
			// Print if:
			//    (1) depth>0 (the default running mode),
			//    (2) depth==0 and the user requested to print zero covered regions (_bedGraphAll)
			if ( (lastDepth != -1) && (lastDepth > 0 || _bedGraphAll) ) {
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
	if ( (lastDepth != -1) && (lastDepth > 0 || _bedGraphAll) ) {
		cout << chrom << "\t" << lastStart << "\t" << chromSize << "\t" << lastDepth << endl;
	}
}
