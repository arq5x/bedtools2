/*
 * ContextMerge.cpp
 *
 *  Created on: Mar 26, 2014
 *      Author: nek3d
 */


#include "ContextMerge.h"

ContextMerge::ContextMerge()
{
	setUseMergedIntervals(true);
	setColumnOpsMethods(true);

	//merge has no default columnOps the way map does, so we'll need to clear those.
	_keyListOps->setColumns("");
	_keyListOps->setOperations("");

}

ContextMerge::~ContextMerge()
{

}


bool ContextMerge::parseCmdArgs(int argc, char **argv, int skipFirstArgs)
{
	_argc = argc;
	_argv = argv;
	_skipFirstArgs = skipFirstArgs;
	if (_argc < 2) {
		setShowHelp(true);
		return false;
	}

	setProgram(_programNames[argv[0]]);

	_argsProcessed.resize(_argc - _skipFirstArgs, false);

	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		else if (strcmp(_argv[_i], "-n") == 0) {
			if (!handle_n()) return false;
		}
		else if (strcmp(_argv[_i], "-nms") == 0) {
			if (!handle_nms()) return false;
		}
		else if (strcmp(_argv[_i], "-scores") == 0) {
			if (!handle_scores()) return false;
		}
		else if (strcmp(_argv[_i], "-delim") == 0) {
			if (!handle_delim()) return false;
		}
		else if (strcmp(_argv[_i], "-d") == 0) {
			if (!handle_d()) return false;
		}
		else if (strcmp(_argv[_i], "-s") == 0) {
			if (!handle_s()) return false;
		}
		else if (strcmp(_argv[_i], "-S") == 0) {
			if (!handle_S()) return false;
		}
	}
	return ContextBase::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextMerge::isValidState()
{
	// Special: The merge program does not have default
	//column operations, so if none were entered, disable column ops.
	if (_keyListOps->getColumns().empty() && _keyListOps->getOperations().empty()) {
		setColumnOpsMethods(false);
		delete _keyListOps;
		_keyListOps = NULL;
	}
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


bool ContextMerge::handle_d() {
    if ((_i+1) < _argc) {
    	if (isNumeric(_argv[_i+1])) {
			int dist = str2chrPos(_argv[_i+1]);
			if (dist >=0 ) {
				_maxDistance = dist;
		    	markUsed(_i - _skipFirstArgs);
		        _i++;
		        markUsed(_i - _skipFirstArgs);
				return true;
			}
    	}
    }
	_errorMsg = "\n***** ERROR: -d option must be followed by an integer value *****";
	return false;
}

bool ContextMerge::handle_n()
{
	//This is the same as telling map "-c any -o count"
	_keyListOps->addColumns("1"); //doesn't really matter which column, but the default column
	//for keyListOps is score, which not every record necessarily has.
	_keyListOps->addOperations("count");
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextMerge::handle_nms()
{
	//This is the same as telling map "-c 4 -o collapse"
	_keyListOps->addColumns("4");
	_keyListOps->addOperations("collapse");
	markUsed(_i - _skipFirstArgs);
    return true;
}


bool ContextMerge::handle_scores()
{
    if ((_i+1) < _argc) {
    	_keyListOps->addColumns("5");
    	_keyListOps->addOperations(_argv[_i+1]);
    	markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
	_errorMsg = "\n***** ERROR: -scores option given, but no operation specified. *****";

    return false;
}

bool ContextMerge::handle_s() {
	_desiredStrand = FileRecordMergeMgr::SAME_STRAND_EITHER;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextMerge::handle_S() {
    if ((_i+1) < _argc) {
    	if (_argv[_i+1][0] == '+') {
			_desiredStrand = FileRecordMergeMgr::SAME_STRAND_FORWARD;
    	} else if (_argv[_i+1][0] == '-') {
			_desiredStrand = FileRecordMergeMgr::SAME_STRAND_REVERSE;
    	}
    	markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
	_errorMsg = "\n***** ERROR: -S option must be followed by + or -. *****";
	return false;
}
