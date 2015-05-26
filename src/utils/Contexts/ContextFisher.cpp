/*
 * ContextFisher.cpp
 *
 */

#include "ContextFisher.h"

ContextFisher::ContextFisher() {
	setSortedInput(true);
	setUseMergedIntervals(false); //by default,
	//intervals are merged in Jaccard, but unmerged in fisher.
}

ContextFisher::~ContextFisher() {
} 
bool ContextFisher::parseCmdArgs(int argc, char **argv, int skipFirstArgs)
{
	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		else if (strcmp(_argv[_i], "-exclude") == 0) {
			if (!handle_exclude()) return false;
		}
		else if (strcmp(_argv[_i], "-m") == 0) {
			if (!handle_m()) return false;
		}

	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextFisher::isValidState()
{
	if (!ContextJaccard::isValidState()) {
		return false;
	}
    if (_genomeFile == NULL){
        _errorMsg = "\nERROR*****: specify -g genome file*****\n";
        return false;
    }
	return true;
}

bool ContextFisher::handle_exclude()
{
    if (_argc <= _i+1) {
        _errorMsg = "\n***** ERROR: -exclude option given, but no file specified. *****";
        return false;
    }

    do {
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    } while (_argc > _i+1 && _argv[_i+1][0] != '-');
    setExcludeFile(string(_argv[_i]));
    return true;
}

bool ContextFisher::handle_m()
{
	setUseMergedIntervals(true);
    markUsed(_i - _skipFirstArgs);
    _i++;
	return true;
}

