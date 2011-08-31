/*****************************************************************************
  fileType.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef FILETYPE_H
#define FILETYPE_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>

using namespace std;

/*****************************************************************************
  Convenience functions to detect whether a given file is
  "regular" and/or "gzipped".

  Kindly contributed by Assaf Gordon.
******************************************************************************/
string string_error(int errnum);
bool isRegularFile(const string& filename);
bool isGzipFile(istream *file);

#endif /* FILETYPE_H */
