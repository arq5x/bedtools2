/*****************************************************************************
fileType.cpp

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include "fileType.h"


/*
returns TRUE if the file is a regular file:
not a pipe/device.

This implies that the file can be opened/closed/seek'd multiple times without losing information
*/
bool isRegularFile(const string& filename) {
    struct stat buf ;
    int i;

    i = stat(filename.c_str(), &buf);
    if (i!=0) {
        cerr << "Error: can't determine file type of '" << filename << "': " << strerror(errno) << endl;
        exit(1);
    }
    if (S_ISREG(buf.st_mode))
        return true;

    return false;
}

/*
returns TRUE if the file has a GZIP header.
Should only be run on regular files.
*/
bool isGzipFile(istream *file) {    
    /*
       11-Sep-2011: 
       We now only peek at the first byte and test for GZIPiness.
       This is because I can only putback() one byte into an istream
       without triggering the "fail" bit.  This was necessary to support
       FIFOs, per version 2.13.0
    */
    struct  {
        unsigned char id1;
//      unsigned char id2;
//      unsigned char cm;
    } gzip_header;

    /*
       26-Dec-2012:
       Just peek at the first byte instead of reading it so that we don't
       affect the istream's failbit.  This modification was wisely proposed
       by John Marshall in response to Issue 30:
       https://github.com/arq5x/bedtools/issues/30
       */
    
    // Check to see if the first byte matches expectations for a
    // GZIP header.
    // Deatils: http://www.gzip.org/zlib/rfc-gzip.html#file-format
    if (file->peek() != 0x1f)
        return false;
    else
        return true;
}
