/*
 * ContextMap.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: arq5x
 */

#include "ContextMap.h"

ContextMap::ContextMap()
{
	// map requires sorted input
	setSortedInput(true);
	setLeftJoin(true);
	setColumnOpsMethods(true);
}

ContextMap::~ContextMap()
{

}


bool ContextMap::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {

	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		if (strcmp(_argv[_i], "-c") == 0) {
			//bypass intersect's use of the -c option, because -c
			//means writeCount for intersect, but means columns for map.
			if (!ContextBase::handle_c()) return false;
		}

	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextMap::isValidState()
{
	if (!ContextIntersect::isValidState()) {
		return false;
	}

	// Multiple databases are currently not supported
	if (getNumDatabaseFiles() > 1) {
		_errorMsg = "\n***** ERROR: multiple database files currently not supported for map. *****";
		return false;
	}
	return true;
}
