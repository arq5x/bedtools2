/*****************************************************************************
genomeCoverage.h

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include "genomeFile.h"

#include "BamReader.h"
#include "BamAncillary.h"
#include "BamAux.h"
using namespace BamTools;

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef map<int, DEPTH, less<int> >          depthMap;
typedef map<string, depthMap, less<string> > chromDepthMap;

typedef map<int, unsigned int, less<int> >   histMap;
typedef map<string, histMap, less<string> >  chromHistMap;

//************************************************
// Class methods and elements
//************************************************
class BedGenomeCoverage {

public:

    // constructor
    BedGenomeCoverage(string bedFile, string genomeFile, bool eachBase, bool startSites,
                      bool bedGraph, bool bedGraphAll, int max, bool bamInput, bool obeySplits,
                      bool filterByStrand, string requestedStrand);

    // destructor
    ~BedGenomeCoverage(void);

private:

    // data (parms)
    string _bedFile;
    string _genomeFile;
    bool   _bamInput;
    bool   _eachBase;
    bool   _startSites;
    bool   _bedGraph;
    bool   _bedGraphAll;
    int    _max;
    bool   _obeySplits;
    bool   _filterByStrand;
    string _requestedStrand;

    BedFile    *_bed;
    GenomeFile *_genome;

    // data for internal processing
    chromDepthMap _chromCov;
    string _currChromName ;
    vector<DEPTH> _currChromCoverage;
    chromHistMap _currChromDepthHist;
    int _currChromSize ;
    set<string> _visitedChromosomes;


    // methods
    void CoverageBed();
    void CoverageBam(string bamFile);
    void ReportChromCoverage(const vector<DEPTH> &, const int &chromSize, const string &chrom, chromHistMap&);
    void ReportGenomeCoverage(chromHistMap &chromDepthHist);
    void ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom);
    void ResetChromCoverage();
    void StartNewChrom (const string& chrom);
    void AddCoverage (int start, int end);
    void AddBlockedCoverage(const vector<BED> &bedBlocks);
    void PrintFinalCoverage();
};
