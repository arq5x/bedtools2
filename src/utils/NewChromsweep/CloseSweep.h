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

#include "ContextClosest.h"

class distanceTuple {
public:
	distanceTuple() : _dist(0), _rec(NULL), _isNeg(false) {}
	distanceTuple(int dist, Record *rec, bool isNeg = false) : _dist(dist), _rec(rec), _isNeg(isNeg) {}
	int _dist;
	Record *_rec;
	bool _isNeg;
};

class DistanceTupleSortAscFunctor {
public:
	bool operator()(const distanceTuple & d1, const distanceTuple & d2) const {
		return (d1._dist < d2._dist ? true : (d1._dist == d2._dist ? d1._rec->lessThan(d2._rec) : false));
	}
};


class RecDistList {
public:
    typedef enum { LEFT, OVERLAP, RIGHT } chromDirType;
	RecDistList(int maxSize);
	~RecDistList();
	bool empty() const { return _empty; }
	void clear();
	int uniqueSize() const { return _currNumIdxs; }
	size_t totalSize() const { return _totalRecs; }
	bool addRec(CHRPOS dist, Record *, chromDirType chromDir);
	bool exists(CHRPOS dist) const {
		CHRPOS dummyVal = 0;
		return find(dist, dummyVal);
	}
	typedef pair<chromDirType, Record *> elemPairType;
	typedef vector<elemPairType> elemsType;
	typedef pair<int, int> indexType;

	int getMaxDist() const { return _empty ? 0 : _distIndex[_currNumIdxs-1].first; }
	typedef int constIterType; //used to be a map iter, trying not to change interface too much.
	constIterType begin() const { return 0; }
	constIterType end() const { return _currNumIdxs; }
	int currDist(constIterType iter) const { return _distIndex[iter].first; }
	size_t currNumElems(constIterType iter) const { return allElems(iter)->size(); }
	elemsType *allElems(constIterType iter) const { return _allRecs[_distIndex[iter].second]; }
	CHRPOS getMaxLeftEndPos() const;

private:

	void insert(CHRPOS dist, Record *, chromDirType chromDir);


	//if true, pos will be the idx the distance is at.
	//if false, pos will be the idx to insert at.
	bool find(CHRPOS dist, CHRPOS &pos) const;


	int _kVal; //max unique allowed
	bool _empty;
	int _currNumIdxs;
	int _totalRecs;

	vector<elemsType *> _allRecs;
	indexType * _distIndex;
};

class CloseSweep : public NewChromSweep {
public:
	CloseSweep(ContextClosest *context);
	~CloseSweep(void);
	bool init();
	const vector<CHRPOS> &getDistances() { return _finalDistances; }

private:
   ContextClosest *_context;
   int _kClosest; // how many closest hits we want to each query.
	vector<RecDistList *> _minUpstreamRecs;
	vector<RecDistList *> _minDownstreamRecs;
	vector<RecDistList *> _overlapRecs;
	vector<CHRPOS> _maxPrevLeftClosestEndPos;
	vector<CHRPOS> _maxPrevLeftClosestEndPosReverse;

	vector<CHRPOS> _finalDistances;


	//
	// Some abbreviations to make the code less miserable.
	//
	bool _sameStrand;
	bool _diffStrand;

	bool _refDist;
	bool _aDist;
	bool _bDist;

	bool _ignoreUpstream;
	bool _ignoreDownstream;

	bool _qForward;
	bool _qReverse;
	bool _dbForward;
	bool _dbReverse;

	ContextClosest::tieModeType _tieMode;
	bool _firstTie;
	bool _lastTie;
	bool _allTies;

	bool allHitsRightOfQueryIgnored(); //true if, no matter what the strands
	// of the hit and query are, we'd ignore the hit so long as it's on the right
	// of the query. Set only during initilization, this is strictly a function
	// of the user provided arguments. Ex: -D ref -id



	//structs to help with finding closest among all of multiple dbs.
	RecordKeyVector _copyRetList;
	vector<CHRPOS> _copyDists;

	//override these methods from chromsweep
	void masterScan(RecordKeyVector &retList);
    void scanCache(int dbIdx, RecordKeyVector &retList);
    bool chromChange(int dbIdx, RecordKeyVector &retList, bool wantScan);


 	typedef enum { IGNORE, DELETE } rateOvlpType;
    rateOvlpType considerRecord(Record *cacheRec, int dbIdx, bool &stopScanning);
    void finalizeSelections(int dbIdx, RecordKeyVector &retList);
    void checkMultiDbs(RecordKeyVector &retList);

    typedef enum { LEFT, OVERLAP, RIGHT } chromDirType;
    typedef enum { UPSTREAM, INTERSECT, DOWNSTREAM } streamDirType;
    typedef enum { NEITHER, FORWARD_ONLY, REVERSE_ONLY, BOTH } purgeDirectionType;

    void setLeftClosestEndPos(int dbIdx);
    bool beforeLeftClosestEndPos(int dbIdx, Record *rec);
    void clearClosestEndPos(int dbIdx);
    int addRecsToRetList(RecDistList::elemsType *recs, CHRPOS currDist, RecordKeyVector &retList);
    void addSingleRec(Record *rec, CHRPOS currDist, int &hitsUsed, RecordKeyVector &retList);
    rateOvlpType tryToAddRecord(Record *cacheRec, CHRPOS dist, int dbIdx, bool &stopScanning, chromDirType chromDir, streamDirType streamDir);
    purgeDirectionType purgePointException();

};


#endif /* CLOSESWEEP_H_ */
