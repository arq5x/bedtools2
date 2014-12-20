/*
 * CloseSweep.h
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

#ifndef CLOSESWEEP_H_
#define CLOSESWEEP_H_

#include "NewChromsweep.h"

class ContextClosest;

class CloseSweep : public NewChromSweep {
public:
	CloseSweep(ContextClosest *context);
	~CloseSweep(void);
	bool init();
	const vector<int> &getDistances() { return _finalDistances; }

private:
   ContextClosest *_context;

	typedef vector<const Record * > distRecVecType;
	vector<distRecVecType *> _minUpstreamRecs;
	vector<int> _minUpstreamDist;
	vector<distRecVecType *> _minDownstreamRecs;
	vector<int> _minDownstreamDist;
	vector<distRecVecType *> _overlapRecs;
	vector<int> _maxPrevLeftClosestEndPos;
	vector<int> _maxPrevLeftClosestEndPosReverse;

	vector<int> _finalDistances;


	//structs to help with finding closest among all of multiple dbs.
	RecordKeyVector _copyRetList;
	vector<int> _copyDists;

	//override these methods from chromsweep
	void masterScan(RecordKeyVector &retList);
    void scanCache(int dbIdx, RecordKeyVector &retList);
    bool chromChange(int dbIdx, RecordKeyVector &retList, bool wantScan);

 	typedef enum { IGNORE, DELETE } rateOvlpType;
    rateOvlpType considerRecord(const Record *cacheRec, int dbIdx, bool &stopScanning);
    void finalizeSelections(int dbIdx, RecordKeyVector &retList);
    void checkMultiDbs(RecordKeyVector &retList);

    void setLeftClosestEndPos(int dbIdx, const Record *rec);
    bool beforeLeftClosestEndPos(int dbIdx, const Record *rec);
    void clearClosestEndPos(int dbIdx);
};


#endif /* CLOSESWEEP_H_ */
