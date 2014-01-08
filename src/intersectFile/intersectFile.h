/*****************************************************************************
  intersectFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef INTERSECTFILE_H
#define INTERSECTFILE_H

using namespace std;

#include "RecordKeyList.h"

using namespace std;

class ContextIntersect;
class BlockMgr;
class RecordOutputMgr;

class FileIntersect {

public:
    FileIntersect(ContextIntersect *context);
    ~FileIntersect(void);

    bool intersectFiles();

private:
    ContextIntersect *_context;
	Record *_queryRec;
	Record *_databaseRec;
	BlockMgr *_blockMgr;
	RecordOutputMgr *_recordOutputMgr;

	void processHits(RecordKeyList &hits);
	bool processSortedFiles();
	bool processUnsortedFiles();

};

#endif /* INTERSECTFILE_H */
