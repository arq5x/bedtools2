/*****************************************************************************
  mapFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef MAPFILE_H
#define MAPFILE_H

#include <sstream>
#include <iomanip>
#include "VectorOps.h"
#include "RecordKeyVector.h"
#include "ContextMap.h"

using namespace std;

class BlockMgr;
class RecordOutputMgr;

class FileMap {

public:
    FileMap(ContextMap *context);
    ~FileMap(void);

    bool mapFiles();

private:
    ContextMap *_context;
    BlockMgr *_blockMgr;
    RecordOutputMgr *_recordOutputMgr;
};

#endif /* MAPFILE_H */
