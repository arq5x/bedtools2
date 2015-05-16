/*
 * mapFile.h
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#ifndef MAPFILE_H_
#define MAPFILE_H_

#include "intersectFile.h"
#include "ContextMap.h"

class MapFile : public IntersectFile {

public:
    MapFile(ContextMap *context);
	virtual bool findNext(RecordKeyVector &hits);
 	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);

protected:
	virtual ContextMap *upCast(ContextBase *context) { return static_cast<ContextMap *>(context); }
};




#endif /* MAPFILE_H_ */
