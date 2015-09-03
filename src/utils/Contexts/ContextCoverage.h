/*
 * ContextCoverage.h
 *
 *  Created on: May 8, 2015
 *      Author: nek3d
 */

#ifndef CONTEXTCOVERAGE_H_
#define CONTEXTCOVERAGE_H_


#include "ContextIntersect.h"

class ContextCoverage : public ContextIntersect {
public:
	ContextCoverage();
	virtual ~ContextCoverage();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
    virtual bool hasIntersectMethods() const { return true; }
    virtual bool isValidState();

    typedef enum { DEFAULT, COUNT, PER_BASE, HIST, MEAN } coverageType;
    coverageType getCoverageType() const { return _coverageType; }

private:
    bool _count;
    bool _perBase;
    bool _showHist;
    bool _mean;
    coverageType _coverageType;

	bool handle_c();
	bool handle_d();
	bool handle_hist();
	bool handle_mean();


};



#endif /* CONTEXTCOVERAGE_H_ */
