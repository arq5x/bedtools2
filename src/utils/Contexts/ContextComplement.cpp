/*
 * ContextComplement.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#include "ContextComplement.h"

ContextComplement::ContextComplement()
{
	setSortedInput(true);
	setUseMergedIntervals(true);
	setExplicitBedOutput(true);

}

ContextComplement::~ContextComplement()
{

}

bool ContextComplement::isValidState()
{
	if (!hasGenomeFile()) {
		_errorMsg = "\n***** ERROR: no -g genome file provided. *****";
		return false;
	}

	return ContextBase::isValidState();
}

bool ContextComplement::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
  for (_i=_skipFirstArgs; _i < argc; _i++) {
    if (isUsed(_i - _skipFirstArgs)) {
      continue;
    }
    else if (strcmp(_argv[_i], "-L") == 0) {
      if (!handle_limit_chroms()) return false;
    }
  }
  return ContextBase::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextComplement::handle_limit_chroms()
{
    _onlyChromsWithBedRecords = true;
    markUsed(_i - _skipFirstArgs);
    return true;
}