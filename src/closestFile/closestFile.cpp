/*
 * newClosestFile.cpp
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

#include "FileRecordMgr.h"
#include "RecordOutputMgr.h"
#include "closestFile.h"
#include "CloseSweep.h"

ClosestFile::ClosestFile(ContextClosest *context)
: _context(context),
 _recordOutputMgr(NULL)
{
	_recordOutputMgr = new RecordOutputMgr();
	_recordOutputMgr->init(_context);
}

ClosestFile::~ClosestFile() {
	delete _recordOutputMgr;
}

bool ClosestFile::getClosest() {
    CloseSweep sweep(_context);
    if (!sweep.init()) {
      return false;
    }
    RecordKeyVector hitSet;
    while (sweep.next(hitSet)) {
    	if (_context->reportDistance()) {
    		_recordOutputMgr->printClosest(hitSet, &(sweep.getDistances()));
    	} else {
    		_recordOutputMgr->printClosest(hitSet, NULL);

    	}
    }
    return true;

}
