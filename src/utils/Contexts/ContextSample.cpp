/*
 * ContextSample.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: nek3d
 */
#include "ContextSample.h"

ContextSample::ContextSample()
{

}

ContextSample::~ContextSample()
{

}

bool ContextSample::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}

        else if (strcmp(_argv[_i], "-s") == 0) {
			if (!handle_s()) return false;
        }
	}
	return ContextBase::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextSample::isValidState()
{
	if (!ContextBase::isValidState()) {
		return false;
	}
	if (_files.size() != 1) {
		_errorMsg = "\n***** ERROR: input file not specified. *****";
		// Allow one and only input file for now
		return false;
	}
	return true;
}


bool ContextSample::handle_s()
{
	_sameStrand = true;
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -s option given, but \"forward\" or \"reverse\" not specified. *****";
		return false;
	}
	if (strcmp(_argv[_i+1], "forward") == 0) {
		_forwardOnly = true;
	} else if (strcmp(_argv[_i+1], "reverse") == 0) {
		_reverseOnly = true;
	} else {
		_errorMsg = "\n***** ERROR: -s option given, but \"forward\" or \"reverse\" not specified. *****";
		return false;
	}
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
    return true;
}


