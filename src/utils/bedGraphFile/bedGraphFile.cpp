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
    if (bedGraphFile == "stdin") {
        _bedGraphStream = &cin;
        return;
    }
    // unzipped, regular
    else if ((isGzipFile(bedGraphFile) == false) && (isRegularFile(bedGraphFile) == true)) {
        _bedGraphStream = new ifstream(bedGraphFile.c_str(), ios::in);

        // open an ifstream
        ifstream bedg(bedGraphFile.c_str(), ios::in);

        // can we open the file?
        if ( !bedg ) {
             cerr << "Error: The requested bedgraph file (" << bedGraphFile << ") could not be opened. Exiting!" << endl;
             exit (1);
         }
         else {
             // if so, close it (this was just a test)
             bedg.close();
             // now set a pointer to the stream so that we
             _bedGraphStream = new ifstream(bedGraphFile.c_str(), ios::in);
         }
     }
     else if ((isGzipFile(bedGraphFile) == true) && (isRegularFile(bedGraphFile) == true)) {

        igzstream bedg(bedGraphFile.c_str(), ios::in);
        if ( !bedg ) {
            cerr << "Error: The requested bedgraph file (" << bedGraphFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            // if so, close it (this was just a test)
            bedg.close();
            // now set a pointer to the stream so that we
            _bedGraphStream = new igzstream(bedGraphFile.c_str(), ios::in);
        }
     }
     else {
         cerr << "Error: Unexpected file type (" << bedGraphFile << "). Exiting!" << endl;
         exit(1);
     }
}


// Close the BEDGRAPH file
void BedGraphFile::Close() {
    if (bedGraphFile != "stdin") {
        if (_bedGraphStream) {
            delete _bedGraphStream;
            _bedGraphStream = NULL ;
        }
    }
}

