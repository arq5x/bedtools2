/*****************************************************************************
genomeCoverage.cpp

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "genomeCoverageBed.h"


BedGenomeCoverage::BedGenomeCoverage(string bedFile, string genomeFile,
                                     bool eachBase, bool startSites, 
                                     bool bedGraph, bool bedGraphAll,
                                     int max, float scale,
                                     bool bamInput, bool obeySplits,
                                     bool filterByStrand, string requestedStrand,
                                     bool only_5p_end, bool only_3p_end,
                                     bool eachBaseZeroBased,
                                     bool add_gb_track_line, string gb_track_line_opts) {

    _bedFile = bedFile;
    _genomeFile = genomeFile;
    _eachBase = eachBase;
    _eachBaseZeroBased = eachBaseZeroBased;
    _startSites = startSites;
    _bedGraph = bedGraph;
    _bedGraphAll = bedGraphAll;
    _max = max;
    _scale = scale;
    _bamInput = bamInput;
    _obeySplits = obeySplits;
    _filterByStrand = filterByStrand;
    _requestedStrand = requestedStrand;
    _only_3p_end = only_3p_end;
    _only_5p_end = only_5p_end;
    _add_gb_track_line = add_gb_track_line;
    _gb_track_line_opts = gb_track_line_opts;
    _currChromName = "";
    _currChromSize = 0 ;

    
    if (_bamInput == false) {
        _genome = new GenomeFile(genomeFile);
    }
    
    PrintTrackDefinitionLine();

    if (_bamInput == false) {
        _bed = new BedFile(bedFile);
        CoverageBed();
    }
    else {
        CoverageBam(_bedFile);
    }
}

void BedGenomeCoverage::PrintTrackDefinitionLine()
{
    //Print Track Definition line (if requested)
    if ( (_bedGraph||_bedGraphAll) && _add_gb_track_line) {
        string line = "track type=bedGraph";
        if (!_gb_track_line_opts.empty()) {
            line += " " ;
            line += _gb_track_line_opts ;
        }
        cout << line << endl;
    }

}


BedGenomeCoverage::~BedGenomeCoverage(void) {
    delete _bed;
    delete _genome;
}


void BedGenomeCoverage::ResetChromCoverage() {
    _currChromName = "";
    _currChromSize = 0 ;
    std::vector<DEPTH>().swap(_currChromCoverage);
}


void BedGenomeCoverage::StartNewChrom(const string& newChrom) {
    // If we've moved beyond the first encountered chromosomes,
    // process the results of the previous chromosome.
    if (_currChromName.length() > 0) {
        ReportChromCoverage(_currChromCoverage, _currChromSize,
                _currChromName, _currChromDepthHist);
    }

    // empty the previous chromosome and reserve new
    std::vector<DEPTH>().swap(_currChromCoverage);

    if (_visitedChromosomes.find(newChrom) != _visitedChromosomes.end()) {
        cerr << "Input error: Chromosome " << _currChromName
             << " found in non-sequential lines. This suggests that the input file is not sorted correctly." << endl;

    }
    _visitedChromosomes.insert(newChrom);

    _currChromName = newChrom;

    // get the current chrom size and allocate space
    _currChromSize = _genome->getChromSize(_currChromName);

    if (_currChromSize >= 0)
        _currChromCoverage.resize(_currChromSize);
    else {
        cerr << "Input error: Chromosome " << _currChromName << " found in your input file but not in your genome file." << endl;
        exit(1);
    }
}


void BedGenomeCoverage::AddCoverage(int start, int end) {
    // process the first line for this chromosome.
    // make sure the coordinates fit within the chrom
    if (start < _currChromSize)
        _currChromCoverage[start].starts++;
    if (end < _currChromSize)
        _currChromCoverage[end].ends++;
    else
        _currChromCoverage[_currChromSize-1].ends++;
}


void BedGenomeCoverage::AddBlockedCoverage(const vector<BED> &bedBlocks) {
    vector<BED>::const_iterator bedItr = bedBlocks.begin();
    vector<BED>::const_iterator bedEnd = bedBlocks.end();
    for (; bedItr != bedEnd; ++bedItr) {
        // the end - 1 must be done because BamAncillary::getBamBlocks
        // returns ends uncorrected for the genomeCoverageBed data structure.
        // ugly, but necessary.
        AddCoverage(bedItr->start, bedItr->end - 1);
    }
}


void BedGenomeCoverage::CoverageBed() {

    BED a;

    ResetChromCoverage();

    _bed->Open();
    while (_bed->GetNextBed(a)) {
        if (_bed->_status == BED_VALID) {
            if (_filterByStrand == true) {
                if (a.strand.empty()) {
                    cerr << "Input error: Interval is missing a strand value on line " << _bed->_lineNum << "." <<endl;
                    exit(1);
                }
                if ( ! (a.strand == "-" || a.strand == "+") ) {
                    cerr << "Input error: Invalid strand value (" << a.strand << ") on line " << _bed->_lineNum << "." << endl;
                    exit(1);
                }
                // skip if the strand is not what the user requested.
                if (a.strand != _requestedStrand)
                    continue;
            }

            // are we on a new chromosome?
            if (a.chrom != _currChromName)
                StartNewChrom(a.chrom);

            if (_obeySplits == true) {
                bedVector bedBlocks; // vec to store the discrete BED "blocks"
                GetBedBlocks(a, bedBlocks);
                AddBlockedCoverage(bedBlocks);
            }
            else if (_only_5p_end) {
                int pos = ( a.strand=="+" ) ? a.start : a.end-1;
                AddCoverage(pos,pos);
            }
            else if (_only_3p_end) {
                int pos = ( a.strand=="-" ) ? a.start : a.end-1;
                AddCoverage(pos,pos);
            }
            else
                AddCoverage(a.start, a.end-1);
        }
    }
    _bed->Close();
    PrintFinalCoverage();
}


void BedGenomeCoverage::PrintFinalCoverage()
{


    // process the results of the last chromosome.
    ReportChromCoverage(_currChromCoverage, _currChromSize,
            _currChromName, _currChromDepthHist);
    if (_eachBase == false && _bedGraph == false && _bedGraphAll == false) {
        ReportGenomeCoverage(_currChromDepthHist);
    }
}


void BedGenomeCoverage::CoverageBam(string bamFile) {

    ResetChromCoverage();

    // open the BAM file
    BamReader reader;
    reader.Open(bamFile);

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // load the BAM header references into a BEDTools "genome file"
    _genome = new GenomeFile(refs);
    // convert each aligned BAM entry to BED
    // and compute coverage on B
    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
        // skip if the read is unaligned
        if (bam.IsMapped() == false)
            continue;

        // skip if we care about strands and the strand isn't what
        // the user wanted
        if ( (_filterByStrand == true) &&
             ((_requestedStrand == "-") != bam.IsReverseStrand()) )
            continue;

        // extract the chrom, start and end from the BAM alignment
        string chrom(refs.at(bam.RefID).RefName);
        CHRPOS start = bam.Position;
        CHRPOS end = bam.GetEndPosition(false, false) - 1;
        
        // are we on a new chromosome?
        if ( chrom != _currChromName )
            StartNewChrom(chrom);

        // add coverage accordingly.
        if (!_only_5p_end && !_only_3p_end) {
            bedVector bedBlocks;
            // we always want to split blocks when a D CIGAR op is found.
            // if the user invokes -split, we want to also split on N ops.
            if (_obeySplits) { // "D" true, "N" true
                GetBamBlocks(bam, refs.at(bam.RefID).RefName, bedBlocks, true, true);
            }
            else { // "D" true, "N" false
                GetBamBlocks(bam, refs.at(bam.RefID).RefName, bedBlocks, true, false);
            }
            AddBlockedCoverage(bedBlocks);
        }
        else if (_only_5p_end) {
            int pos = ( !bam.IsReverseStrand() ) ? start : end;
            AddCoverage(pos,pos);
        }
        else if (_only_3p_end) {
            int pos = ( bam.IsReverseStrand() ) ? start : end;
            AddCoverage(pos,pos);
        }
    }
    // close the BAM
    reader.Close();
    PrintFinalCoverage();
}


void BedGenomeCoverage::ReportChromCoverage(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom, chromHistMap &chromDepthHist) {

    if (_eachBase) {
        int depth = 0; // initialize the depth
        int offset = (_eachBaseZeroBased)?0:1;
        for (int pos = 0; pos < chromSize; pos++) {

            depth += chromCov[pos].starts;
            // report the depth for this position.
            if (depth>0 || !_eachBaseZeroBased)
                cout << chrom << "\t" << pos+offset << "\t" << depth * _scale << endl;
            depth = depth - chromCov[pos].ends;
        }
    }
    else if (_bedGraph == true || _bedGraphAll == true) {
        ReportChromCoverageBedGraph(chromCov, chromSize, chrom);
    }
    else {

        int depth = 0; // initialize the depth

        for (int pos = 0; pos < chromSize; pos++) {

            depth += chromCov[pos].starts;

            // add the depth at this position to the depth histogram
            // for this chromosome. if the depth is greater than the
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
        histMap::const_iterator depthIt = chromDepthHist[chrom].begin();
        histMap::const_iterator depthEnd = chromDepthHist[chrom].end();
        for (; depthIt != depthEnd; ++depthIt) {
            int depth = depthIt->first;
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
        genomeSize += _genome->getChromSize(chrom);
        // if there were no reads for a give chromosome, then
        // add the length of the chrom to the 0 bin.
        if ( chromDepthHist.find(chrom) == chromDepthHist.end() ) {
            chromDepthHist[chrom][0] += _genome->getChromSize(chrom);
        }
    }

    histMap genomeHist; // depth histogram for the entire genome

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


void BedGenomeCoverage::ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom) {

    int depth = 0; // initialize the depth
    int lastStart = -1;
    int lastDepth = -1;

    for (int pos = 0; pos < chromSize; pos++) {
        depth += chromCov[pos].starts;

        if (depth != lastDepth) {
            // Coverage depth has changed, print the last interval coverage (if any)
            // Print if:
            // (1) depth>0 (the default running mode),
            // (2) depth==0 and the user requested to print zero covered regions (_bedGraphAll)
            if ( (lastDepth != -1) && (lastDepth > 0 || _bedGraphAll) ) {
                cout << chrom << "\t" << lastStart << "\t" << pos << "\t" << lastDepth * _scale << endl;
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
        cout << chrom << "\t" << lastStart << "\t" << chromSize << "\t" << lastDepth * _scale << endl;
    }
}
