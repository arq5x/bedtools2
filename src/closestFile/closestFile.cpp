/*
 * closestFile.cpp
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#include "closestFile.h"
#include "CloseSweep.h"

ClosestFile::ClosestFile(ContextClosest *context)
: IntersectFile(context)
{

}

bool ClosestFile::findNext(RecordKeyVector &hits)
{
	return nextSortedFind(hits);
}

void ClosestFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	if (upCast(_context)->reportDistance()) {
		outputMgr->printClosest(hits, &(upCastSweep()->getDistances()));
	} else {
		outputMgr->printClosest(hits, NULL);
	}
}



void ClosestFile::makeSweep() {
	_sweep = new CloseSweep(upCast(_context));
}
