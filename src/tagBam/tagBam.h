/*****************************************************************************
  tagBam.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef TAGBAM_H
#define TAGBAM_H

#include "bedFile.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "BamAncillary.h"
using namespace BamTools;

#include "bedFile.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class TagBam {

public:

    // constructor
    TagBam(const string &bamFile, const vector<string> &annoFileNames,
                const vector<string> &annoLabels, const string &tag, 
                bool useNames, bool useScores, bool useIntervals, bool sameStrand, 
                bool diffStrand, float overlapFraction);

    // destructor
    ~TagBam(void);

    // annotate the BAM file with all of the annotation files.
    void Tag();

private:

    // input files.
    string _bamFile;
    vector<string> _annoFileNames;
    vector<string> _annoLabels;
        
    string _tag;

    // instance of a bed file class.
    vector<BedFile*> _annoFiles;

    // should we use the name field from the annotation files?
    bool _useNames;
    bool _useScores;
    bool _useIntervals;
    
    // do we care about strandedness when tagging?
    bool _sameStrand;
    bool _diffStrand;
    float _overlapFraction;

    // private function for reporting coverage information
    void ReportAnnotations();

    void OpenAnnoFiles();

    void CloseAnnoFiles();

};
#endif /* TAGBAM_H */
