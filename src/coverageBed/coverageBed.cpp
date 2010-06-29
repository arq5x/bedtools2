/*****************************************************************************
  coverageBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "coverageBed.h"

// build
BedCoverage::BedCoverage(string &bedAFile, string &bedBFile, bool &forceStrand, 
                         bool &writeHistogram, bool &bamInput, bool &obeySplits) {
	
	_bedAFile       = bedAFile;
	_bedBFile       = bedBFile;
	
	_bedA           = new BedFile(bedAFile);
	_bedB           = new BedFile(bedBFile);
	
	_forceStrand    = forceStrand;
    _obeySplits     = obeySplits;
	_writeHistogram = writeHistogram;
	_bamInput       = bamInput;
	
	if (_bedA->bedFile != "stdin") {   // process a file
		if (_bamInput == false) { //bed/gff
			CollectCoverageBed();
		}
		else {
			CollectCoverageBam(_bedA->bedFile);
		}
	}
	else {   // process stdin
		if (_bamInput == false) 
			CollectCoverageBed();
		else {
			CollectCoverageBam("stdin");	
		}
	}
}

// destroy
BedCoverage::~BedCoverage(void) {
	delete _bedA;
	delete _bedB;
}


void BedCoverage::CollectCoverageBed() {
	
	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	_bedB->loadBedCovFileIntoMap();

	int lineNum = 0;					// current input line number
	BED a, nullBed;	
	BedLineStatus bedStatus;
	
	_bedA->Open();	
	// process each entry in A
	while ((bedStatus = _bedA->GetNextBed(a, lineNum)) != BED_INVALID) {
		if (bedStatus == BED_VALID) {
		    // process the BED entry as a single block
            if (_obeySplits == false)
    			_bedB->countHits(a, _forceStrand);
			// split the BED into discrete blocksand process each independently.
			else {
			    bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
                splitBedIntoBlocks(a, lineNum, bedBlocks);
                
                vector<BED>::const_iterator bedItr  = bedBlocks.begin();
            	vector<BED>::const_iterator bedEnd  = bedBlocks.end();
            	for (; bedItr != bedEnd; ++bedItr)
        	        _bedB->countHits(*bedItr, _forceStrand);
			}
			a = nullBed;
		}
	}	
	_bedA->Close();
	
	// report the coverage (summary or histogram) for BED B.
	ReportCoverage();					
}


void BedCoverage::CollectCoverageBam(string bamFile) {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	_bedB->loadBedCovFileIntoMap();

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
		    // treat the BAM alignment as a single "block"
		    if (_obeySplits == false) {
		        // construct a new BED entry from the current BAM alignment.	
    			BED a;
    			a.chrom  = refs.at(bam.RefID).RefName;
    			a.start  = bam.Position;
    			a.end    = bam.GetEndPosition(false);
    			a.strand = "+"; 
    			if (bam.IsReverseStrand()) a.strand = "-";
    			
    			_bedB->countHits(a, _forceStrand);
			}
			// split the BAM alignment into discrete blocks and
			// look for overlaps only within each block.
			else {
			    // vec to store the discrete BED "blocks" from a
			    bedVector bedBlocks;
                // since we are counting coverage, we do want to split blocks when a 
                // deletion (D) CIGAR op is encountered (hence the true for the last parm)
                getBamBlocks(bam, refs, bedBlocks, true);   
                
                vector<BED>::const_iterator bedItr  = bedBlocks.begin();
            	vector<BED>::const_iterator bedEnd  = bedBlocks.end();
            	for (; bedItr != bedEnd; ++bedItr) {
            	    _bedB->countHits(*bedItr, _forceStrand);
        	    }
			} 	
		}
	}
	// report the coverage (summary or histogram) for BED B.
	ReportCoverage();
	// close the BAM file	
	reader.Close();
}


void BedCoverage::ReportCoverage() {

	map<unsigned int, unsigned int> allDepthHist;
	unsigned int totalLength = 0;

	masterBedCovMap::const_iterator chromItr = _bedB->bedCovMap.begin();
	masterBedCovMap::const_iterator chromEnd = _bedB->bedCovMap.end();
	for (; chromItr != chromEnd; ++chromItr) {
	
		binsToBedCovs::const_iterator binItr = chromItr->second.begin();
		binsToBedCovs::const_iterator binEnd = chromItr->second.end();
		for (; binItr != binEnd; ++binItr) {

			vector<BEDCOV>::const_iterator bedItr = binItr->second.begin();
			vector<BEDCOV>::const_iterator bedEnd = binItr->second.end();
			for (; bedItr != bedEnd; ++bedItr) {
									
				int zeroDepthCount = 0;
				int depth          = 0;
				int start          = min(bedItr->minOverlapStart, bedItr->start);
				
				// track the numnber of bases in the feature covered by
				// 0, 1, 2, ... n features in A
				map<unsigned int, unsigned int> depthHist;
				map<unsigned int, DEPTH>::const_iterator depthItr;
				
				for (CHRPOS pos = start+1; pos <= bedItr->end; pos++) {
					
					depthItr = bedItr->depthMap.find(pos);
					
					if (depthItr != bedItr->depthMap.end()) {
						depth += depthItr->second.starts;
						if ((pos > bedItr->start) && (pos <= bedItr->end)) {	
							if (depth == 0) zeroDepthCount++;
							depthHist[depth]++;
							allDepthHist[depth]++;
						}
						depth = depth - depthItr->second.ends;
					}
					else {
						if ((pos > bedItr->start) && (pos <= bedItr->end)) {	
							if (depth == 0) zeroDepthCount++;
							depthHist[depth]++;
							allDepthHist[depth]++;
						}
					}
				}

				// Report the coverage for the current interval.
				CHRPOS length   = bedItr->end - bedItr->start;
				totalLength += length;
				
				int nonZeroBases   = (length - zeroDepthCount);
				float fractCovered = (float) nonZeroBases / length;
				
				if (_writeHistogram == false) {
					_bedB->reportBedTab(*bedItr);
					printf("%d\t%d\t%d\t%0.7f\n", bedItr->count, nonZeroBases, length, fractCovered);
				}
				else {
					map<unsigned int, unsigned int>::const_iterator histItr = depthHist.begin();
					map<unsigned int, unsigned int>::const_iterator histEnd = depthHist.end();
					for (; histItr != histEnd; ++histItr) {
						float fractAtThisDepth = (float) histItr->second / length;
						_bedB->reportBedTab(*bedItr);
						printf("%d\t%d\t%d\t%0.7f\n", histItr->first, histItr->second, length, fractAtThisDepth);
					}
				}
			}
		}
	}
	if (_writeHistogram == true) {
		map<unsigned int, unsigned int>::const_iterator histItr = allDepthHist.begin();
		map<unsigned int, unsigned int>::const_iterator histEnd = allDepthHist.end();
		for (; histItr != histEnd; ++histItr) {
			float fractAtThisDepth = (float) histItr->second / totalLength;
			printf("all\t%d\t%d\t%d\t%0.7f\n", histItr->first, histItr->second, totalLength, fractAtThisDepth);
		}
	}
}


