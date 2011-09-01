/*****************************************************************************
  bedGraphFile.cpp

  (c) 2010 - Assaf Gordon
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedGraphFile.h"
#include <sstream>

// Constructor
BedGraphFile::BedGraphFile(string &_file) :
    bedGraphFile(_file),
    _bedGraphStream(NULL)
{}


// Destructor
BedGraphFile::~BedGraphFile() {
    Close();
}


// Open the BEDGRAPH file
void BedGraphFile::Open() {
    if (bedGraphFile == "stdin" || bedGraphFile == "-") {
        _bedGraphStream = &cin;
    }
    else {
        _bedGraphStream = new ifstream(bedGraphFile.c_str(), ios::in);

        if (isGzipFile(_bedGraphStream) == true) {
            delete _bedGraphStream;
            _bedGraphStream = new igzstream(bedGraphFile.c_str(), ios::in);
        }
        // can we open the file?
        if ( !(_bedGraphStream->good()) ) {
            cerr << "Error: The requested bed file (" << bedGraphFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
    }
}


// Close the BEDGRAPH file
void BedGraphFile::Close() {
    if (bedGraphFile != "stdin" && bedGraphFile != "-") {
        if (_bedGraphStream) {
            delete _bedGraphStream;
            _bedGraphStream = NULL ;
        }
    }
}

