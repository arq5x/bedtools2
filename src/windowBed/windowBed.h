/*****************************************************************************
  windowBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef WINDOWBED_H
#define WINDOWBED_H

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
using namespace BamTools;

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedWindow {

public:

    // constructor
    BedWindow(string bedAFile, string bedBFile, int leftSlop, int rightSlop,
              bool anyHit, bool noHit, bool writeCount, bool strandWindows,
              bool matchOnSameStrand, bool matchOnDiffStrand, bool bamInput,
              bool bamOutput, bool isUncompressedBam, bool printHeader);

    // destructor
    ~BedWindow(void);

private:

    string _bedAFile;
    string _bedBFile;
    bool _anyHit;
    bool _writeCount;
    int _leftSlop;
    int _rightSlop;
    bool _noHit;
    bool _strandWindows;
    bool _matchOnSameStrand;
    bool _matchOnDiffStrand;
    bool _bamInput;
    bool _bamOutput;
    bool _isUncompressedBam;
    bool _printHeader;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    // methods
    void WindowIntersectBed();
    void WindowIntersectBam(string bamFile);
    void FindWindowOverlaps(const BED &a, vector<BED> &hits);
    bool FindOneOrMoreWindowOverlaps(const BED &a);
    void AddWindow(const BED &a, CHRPOS &fudgeStart, CHRPOS &fudgeEnd);

};
#endif /* WINDOWBED_H */
