/*****************************************************************************
  tabFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef TABFILE_H
#define TABFILE_H

#include "gzstream.h"
#include "fileType.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

// enum to flag the state of a given line in a TAB file.
enum TabLineStatus
{
    TAB_INVALID = -1,
    TAB_HEADER  = 0,
    TAB_BLANK   = 1,
    TAB_VALID   = 2
};

typedef vector<string> TAB_FIELDS;

//************************************************
// TabFile Class methods and elements
//************************************************
class TabFile {

public:

    // Constructor
    TabFile(const string &tabFile);

    // Destructor
    ~TabFile(void);

    // Open a TAB file for reading (creates an istream pointer)
    void Open(void);

    // Close an opened TAB file.
    void Close(void);

    // Get the next TAB entry in an opened TAB file.
    TabLineStatus GetNextTabLine (TAB_FIELDS &tab, int &lineNum);

private:

    // data
    istream *_tabStream;
    string _tabFile;

    // methods
    inline TabLineStatus parseTabLine (const vector<string> &lineVector, int &lineNum) {
        // bail out if we have a blank line
        if (lineVector.size() == 0)
            return TAB_BLANK;
        // real line with data
        if (lineVector[0][0] != '#') {
            return TAB_VALID;
        }
        // comment or header line
        else {
            lineNum--;
            return TAB_HEADER;
        }
        // default
        return TAB_INVALID;
    }
};

#endif /* TABFILE_H */
