/*****************************************************************************
  mergeFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "mergeFile.h"


MergeFile::MergeFile(ContextMerge *context)
: _context(context),
  _recordOutputMgr(NULL)
{
	_recordOutputMgr = new RecordOutputMgr();
	_recordOutputMgr->init(_context);
}

MergeFile::~MergeFile()
{
	delete _recordOutputMgr;
	_recordOutputMgr = NULL;
}

bool MergeFile::merge()
{
    RecordKeyList hitSet;
    FileRecordMgr *frm = _context->getFile(0);
    while (!frm->eof()) {
    	Record *key = frm->getNextRecord(&hitSet);
    	if (key == NULL) continue;
		_recordOutputMgr->printRecord(hitSet.getKey(), _context->getColumnOpsVal(hitSet));
    }
    return true;
}
