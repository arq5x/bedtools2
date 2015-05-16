/*
 * mapFile.cpp
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#include "mapFile.h"

MapFile::MapFile(ContextMap *context)
: IntersectFile(context)
{

}

bool MapFile::findNext(RecordKeyVector &hits)
{
	if (nextSortedFind(hits)) {
		 checkSplits(hits);
		 return true;
	}
	return false;
}

void MapFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	outputMgr->printRecord(hits.getKey(), _context->getColumnOpsVal(hits));
}

