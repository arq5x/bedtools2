/*
 * ContextQc.cpp
 *
 *  Created on: Jan 15, 2018
 *      Author: Aaron Quinlan
 */
#include "ContextQc.h"

ContextQc::ContextQc()
{
}

ContextQc::~ContextQc()
{

}

bool ContextQc::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {

	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
	}
	return ContextBase::parseCmdArgs(argc, argv, 1);
}


bool ContextQc::isValidState()
{
	// if (!hasGenomeFile()) {
	// 	_errorMsg = "\n***** ERROR: no -g genome file provided. *****";
	// 	return false;
	// }
	if (!ContextBase::isValidState()) {
		return false;
	}
	return true;
}