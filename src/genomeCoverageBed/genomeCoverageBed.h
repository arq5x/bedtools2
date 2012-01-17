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

#include "BlockedIntervals.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
using namespace BamTools;

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef map<int, DEPTH, less<int> > depthMap;
typedef map<string, depthMap, less<string> > chromDepthMap;

typedef map<int, unsigned int, less<int> > histMap;
typedef map<string, histMap, less<string> > chromHistMap;

//************************************************
// Class methods and elements
//************************************************
class BedGenomeCoverage {

public:

    // constructor
    BedGenomeCoverage(string bedFile, string genomeFile, 
                      bool eachBase, bool startSites,
                      bool bedGraph, bool bedGraphAll, 
                      int max, float scale, 
                      bool bamInput, bool obeySplits,
                      bool filterByStrand, string requestedStrand,
                      bool only_5p_end, bool only_3p_end,
                      bool eachBaseZeroBased,
                      bool add_gb_track_line, string gb_track_line_opts);

    // destructor
    ~BedGenomeCoverage(void);

private:

    // data (parms)
    string _bedFile;
    string _genomeFile;
    bool _bamInput;
    bool _eachBase;
    bool _eachBaseZeroBased;
    bool _startSites;
    bool _bedGraph;
    bool _bedGraphAll;
    int _max;
    float _scale;
    bool _obeySplits;
    bool _filterByStrand;
    bool _only_5p_end;
    bool _only_3p_end;
    bool _add_gb_track_line;
    string _gb_track_line_opts;
    string _requestedStrand;

    BedFile *_bed;
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
    void LoadBamHeaderIntoGenomeFile(const string &bamFile);
    void ReportChromCoverage(const vector<DEPTH> &, const int &chromSize, const string &chrom, chromHistMap&);
    void ReportGenomeCoverage(chromHistMap &chromDepthHist);
    void ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom);
    void ResetChromCoverage();
    void StartNewChrom (const string& chrom);
    void AddCoverage (int start, int end);
    void AddBlockedCoverage(const vector<BED> &bedBlocks);
    void PrintFinalCoverage();
    void PrintTrackDefinitionLine();
};

