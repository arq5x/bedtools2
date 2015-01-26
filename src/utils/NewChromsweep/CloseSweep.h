/*
 * CloseSweep.h
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

#ifndef CLOSESWEEP_H_
#define CLOSESWEEP_H_

#include "NewChromsweep.h"
#include <list>
#include <set>

class ContextClosest;

class distanceTuple {
public:
	distanceTuple() : _dist(0), _rec(NULL), _isNeg(false) {}
	distanceTuple(int dist, const Record *rec, bool isNeg = false) : _dist(dist), _rec(rec), _isNeg(isNeg) {}
//	bool operator < (const distanceTuple & other) const { return (_dist < other._dist); }
	int _dist;
	const Record *_rec;
	bool _isNeg;
};

class DistanceTupleSortAscFunctor {
public:
	bool operator()(const distanceTuple & d1, const distanceTuple & d2) const {
		return ((d1._dist < d2._dist) ? true : (d1._dist == d2._dist ? (d1._isNeg && !d2._isNeg) : false)) ; }
};


class RecDistList {
public:
    typedef enum { LEFT, OVERLAP, RIGHT } chromDirType;
	RecDistList(int maxSize) : _maxUniqueAllowed(maxSize) {}
	~RecDistList() { clear(); }
	bool empty() const { return _recs.empty(); }
	void clear();
	int uniqueSize() const { return _recs.size(); }
	size_t totalSize() const { return _totalRecs; }
	bool addRec(int dist, const Record *, chromDirType chromDir);
	bool exists(int dist) const { return (_recs.find(dist) != _recs.end()); }
	int furtherestDistance() const { return _recs.rbegin()->first; }

	typedef vector<pair<chromDirType, const Record *> >elemsType;
	typedef map<int, elemsType *> distRecsType;
	typedef distRecsType::iterator iterType;
	typedef distRecsType::reverse_iterator revIterType;
	typedef distRecsType::const_iterator constIterType;
	typedef distRecsType::const_reverse_iterator constRevIterType;

	int getMaxDist() const { return _recs.empty() ? 0 : _recs.rbegin()->first; }
	constIterType begin() const { return _recs.begin(); }
	constIterType end() const { return _recs.end(); }
	int currDist(constIterType iter) const { return iter->first; }
	size_t currNumElems(constIterType iter) const { return iter->second->size(); }
	const Record *firstElem(constIterType iter) const { return (iter->second->at(0)).second; }
	const Record *lastElem(constIterType iter) const { return (iter->second->at(iter->second->size()-1)).second; }
	const elemsType *allElems(constIterType iter) const { return iter->second; }
	int getMaxLeftEndPos() const;

private:
	void insert(int dist, const Record *, chromDirType chromDir);
	distRecsType _recs;
	int _maxUniqueAllowed;
	int _totalRecs;

};

class CloseSweep : public NewChromSweep {
public:
	CloseSweep(ContextClosest *context);
	~CloseSweep(void);
	bool init();
	const vector<int> &getDistances() { return _finalDistances; }

private:
   ContextClosest *_context;
   int _kClosest; // how many closest hits we want to each query.
	vector<RecDistList *> _minUpstreamRecs;
	vector<RecDistList *> _minDownstreamRecs;
	vector<RecDistList *> _overlapRecs;
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

    typedef enum { LEFT, OVERLAP, RIGHT } chromDirType;
    typedef enum { UPSTREAM, INTERSECT, DOWNSTREAM } streamDirType;

    void setLeftClosestEndPos(int dbIdx);
    bool beforeLeftClosestEndPos(int dbIdx, const Record *rec);
    void clearClosestEndPos(int dbIdx);
    bool canStopScan(const Record *cacheRec, bool ignored, streamDirType streamDir);
    int addRecsToRetList(const RecDistList::elemsType *recs, int currDist, RecordKeyVector &retList);
    void addSingleRec(const Record *rec, int currDist, int &hitsUsed, RecordKeyVector &retList);
    rateOvlpType tryToAddRecord(const Record *cacheRec, int dist, int dbIdx, bool &stopScanning, chromDirType chromDir, streamDirType streamDir);
    bool purgePointException(int dbIdx);

};


#endif /* CLOSESWEEP_H_ */
