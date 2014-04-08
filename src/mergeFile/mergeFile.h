/*****************************************************************************
  mergeFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef MERGE_FILE_H_
#define MERGE_FILE_H_

//************************************************
// Class methods and elements
//************************************************

#include "ContextMerge.h"
#include "RecordOutputMgr.h"

class MergeFile {

public:
  MergeFile(ContextMerge *context);
  ~MergeFile();

  bool merge();

private:
  ContextMerge *_context;
  RecordOutputMgr *_recordOutputMgr;
    
};

#endif
