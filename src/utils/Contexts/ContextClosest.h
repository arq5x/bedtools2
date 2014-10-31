/*
 * ContextClosest.h
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */


#ifndef CONTEXTCLOSEST_H_
#define CONTEXTCLOSEST_H_

#include "ContextIntersect.h"

class ContextClosest : public ContextIntersect {
public:
	ContextClosest();
	virtual ~ContextClosest();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
    virtual bool hasIntersectMethods() const { return true; }
    virtual bool isValidState();

    bool hasTieMode() const { return _haveTieMode; }
    bool ignoreOverlaps() const { return _ignoreOverlaps; }
    bool ignoreUpstream() const { return _ignoreUpstream; }
    bool ignoreDownstream() const { return _ignoreDownstream; }
    bool reportDistance() const { return _reportDistance; }
    bool signDistance() const { return _signDistance; }
    bool hasStrandedDistMode() const { return _haveStrandedDistMode; }
    bool diffNames() const { return _diffNames; }

    typedef enum { FIRST_TIE, LAST_TIE, ALL_TIES} tieModeType;
    tieModeType getTieMode() const { return _tieMode; }

    typedef enum { REF_DIST, A_DIST, B_DIST} strandedDistanceModeType;
    strandedDistanceModeType getStrandedDistMode() const { return _strandedDistMode; }

    typedef enum { EACH_DB, ALL_DBS } multiDbModeType;
    multiDbModeType getMultiDbMode() const { return _multiDbMode; }

private:
    bool _haveTieMode;
    bool _ignoreOverlaps;
    bool _ignoreUpstream;
    bool _ignoreDownstream;
    bool _reportDistance;
    bool _signDistance;
    bool _haveStrandedDistMode;
    bool _diffNames;
    tieModeType _tieMode;
    strandedDistanceModeType _strandedDistMode;
    multiDbModeType _multiDbMode;

    bool handle_d();
    bool handle_D();
    bool handle_io();
    bool handle_iu();
    bool handle_id();
    bool handle_N();
    bool handle_t();
    bool handle_mdb();

};


#endif /* CONTEXTCLOSEST_H_ */
