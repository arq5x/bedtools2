/*****************************************************************************
  mapBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "mapBed.h"

const int PRECISION = 21;
double GetUserColumn(const string s);

// Constructor
BedMap::BedMap(string bedAFile, string bedBFile, int column, string operation,
               float overlapFraction, bool sameStrand, 
               bool diffStrand, bool reciprocal, 
               bool printHeader) {

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _column              = column - 1;
    _operation           = operation;
    _overlapFraction     = overlapFraction;
    _sameStrand          = sameStrand;
    _diffStrand          = diffStrand;
    _reciprocal          = reciprocal;
    _printHeader         = printHeader;
    
    Map();
}

// Destructor
BedMap::~BedMap(void) {
}

void BedMap::Map() {

    // create new BED file objects for A and B
    _bedA = new BedFile(_bedAFile);
    _bedB = new BedFile(_bedBFile);

    // use the chromsweep algorithm to detect overlaps on the fly.
    ChromSweep sweep = ChromSweep(_bedB, _bedA, _sameStrand, _diffStrand, _printHeader);

    pair<BED, vector<BED> > hit_set;
    hit_set.second.reserve(100000);
    while (sweep.Next(hit_set)) {
        string result = ApplyHits(hit_set.first, hit_set.second);
        _bedB->reportBedTab(hit_set.first);
        printf("%s\n", result.c_str());
    }
}

string BedMap::MapHits(const BED &a, const vector<BED> &hits) {

    vector<string> data;
    vector<double> dataF;
    data.reserve(hits.size());
    dataF.reserve(hits.size());
    for (size_t i = 0; i < hits.size(); ++i) {
        try {
            data.push_back(hits[i].fields.at(_column));
        }
        catch(std::out_of_range& e) {
            cerr << endl << "*****" << endl 
                 << "*****ERROR: requested column ("
                 << _column + 1
                 << ") exceeds the number of columns in file at line "
                 << _bedA->_lineNum << ". Exiting." 
                 << endl << "*****" << endl;
            exit(1);
        }
    }
    transform(data.begin(), data.end(), back_inserter(dataF), GetUserColumn);

    // sum
    double total = accumulate(dataF.begin(), dataF.end(), 0.0);
    ostringstream output;
    output << setprecision (PRECISION) << total;
    return output.str();
}


double GetUserColumn(const string s) {
    std::istringstream i(s);
    double x;
    if (!(i >> x)) {
        cerr << "Error: Could not properly convert string to numeric (\"" + i.str() + "\")" << endl;
        exit(1);
    }
    return x;
}


