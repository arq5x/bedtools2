/*****************************************************************************
  tabFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "tabFile.h"

/*******************************************
Class methods
*******************************************/

// Constructor
TabFile::TabFile(const string &tabFile)
: _tabFile(tabFile)
{}

// Destructor
TabFile::~TabFile(void) {
}

void TabFile::Open(void) {
    if (_tabFile == "stdin" || _tabFile == "-") {
        _tabStream = &cin;
    }
    else {
        _tabStream = new ifstream(_tabFile.c_str(), ios::in);
        
        if( isGzipFile(_tabStream) ) {
            delete _tabStream;
            _tabStream = new igzstream(_tabFile.c_str(), ios::in);
        }
        if ( _tabStream->fail() ) {
            cerr << "Error: The requested file (" 
                 << _tabFile
                 << ") " 
                 << "could not be opened. "
                 << "Error message: ("
                 << strerror(errno)
                 << "). Exiting!" << endl;
            exit (1);
        }
    }
}


// Close the TAB file
void TabFile::Close(void) {
    if (_tabFile != "stdin" && _tabFile != "-") delete _tabStream;
}


TabLineStatus TabFile::GetNextTabLine(TAB_FIELDS &tabFields, int &lineNum) {

    // make sure there are still lines to process.
    // if so, tokenize, return the TAB_FIELDS.
    if (_tabStream->good() == true) {
        string tabLine;
        tabFields.reserve(20);

        // parse the tabStream pointer
        getline(*_tabStream, tabLine);
        
        if (tabLine[tabLine.size()-1] == '\r') {
            tabLine.resize(tabLine.size()-1);
        }
        
        lineNum++;

        // split into a string vector.
        Tokenize(tabLine, tabFields);

        // parse the line and validate it
        return parseTabLine(tabFields, lineNum);
    }

    // default if file is closed or EOF
    return TAB_INVALID;
}
