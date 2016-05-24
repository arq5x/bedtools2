/*
 * ContextCoverage.cpp
 *
 *  Created on: May 8, 2015
 *      Author: nek3d
 */

#include "ContextCoverage.h"

ContextCoverage::ContextCoverage()
: _count(false),
  _perBase(false),
  _showHist(false),
  _coverageType(DEFAULT)
{
	setExplicitBedOutput(true); //do not allow BAM output
	setRunToQueryEnd(true);
}

ContextCoverage::~ContextCoverage() {
}

bool ContextCoverage::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {

	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		if (strcmp(_argv[_i], "-d") == 0) {
			if (!handle_d()) return false;
		}
		if (strcmp(_argv[_i], "-counts") == 0) {
			if (!handle_c()) return false;
		}
		if (strcmp(_argv[_i], "-hist") == 0) {
			if (!handle_hist()) return false;
		}
		if (strcmp(_argv[_i], "-mean") == 0) {
			if (!handle_mean()) return false;
		}

	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextCoverage::isValidState()
{
	if (!ContextIntersect::isValidState()) {
		return false;
	}
	//Can only use one output option. Were two or more set?
	if (((int)_count + (int)_perBase + (int)_showHist) + (int)_mean > 1) {
		_errorMsg = "\n***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****";
		return false;
	}

	return true;
}

bool ContextCoverage::handle_c()
{
	_count = true;
	_coverageType = COUNT;
	return ContextIntersect::handle_c();
}

bool ContextCoverage::handle_d()
{
	_perBase = true;
	_coverageType = PER_BASE;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextCoverage::handle_hist()
{
	_showHist = true;
	_coverageType = HIST;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextCoverage::handle_mean()
{
	_mean = true;
	_coverageType = MEAN;
	markUsed(_i - _skipFirstArgs);
	return true;
}
