/*****************************************************************************
  annotateBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef ANNOTATEBED_H
#define ANNOTATEBED_H

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
class BedAnnotate {

public:

    // constructor
    BedAnnotate(const string &mainFile, const vector<string> &annoFileNames,
                const vector<string> &annoTitles, bool sameStrand, bool diffStrand, bool reportCounts, bool reportBoth);

    // destructor
    ~BedAnnotate(void);

    // annotate the master file with all of the annotation files.
    void AnnotateBed();

private:

    // input files.
    string _mainFile;
    vector<string> _annoFileNames;
    vector<string> _annoTitles;

    // instance of a bed file class.
    BedFile *_bed;
    vector<BedFile*> _annoFiles;

    // do we care about strandedness when counting coverage?
    bool _sameStrand;
    bool _diffStrand;
    
    bool _reportCounts;
    bool _reportBoth;

    // private function for reporting coverage information
    void ReportAnnotations();

    void OpenAnnoFiles();

    void CloseAnnoFiles();

    void PrintHeader();

    void InitializeMainFile();
};
#endif /* ANNOTATEBED_H */
