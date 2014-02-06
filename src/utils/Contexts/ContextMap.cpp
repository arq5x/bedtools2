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

	// default to BED score column
	setColumn(5);
	// default to "sum"
	setColumnOperation("sum");
	// default to "." as a NULL value
	setNullValue('.');
}

ContextMap::~ContextMap()
{

}


bool ContextMap::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
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
        else if (strcmp(_argv[_i], "-o") == 0) {
			if (!handle_o()) return false;
        }
        else if (strcmp(_argv[_i], "-c") == 0) {
			if (!handle_c()) return false;
        }
        else if (strcmp(_argv[_i], "-null") == 0) {
			if (!handle_null()) return false;
        }
	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}


bool ContextMap::isValidState()
{
	if (!ContextIntersect::isValidState()) {
		return false;
	}

    if (getDatabaseFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
         //throw Error
        cerr << endl << "*****" 
             << endl 
             << "***** ERROR: BAM database file not currently supported for the map tool." 
             << endl;
        exit(1);
    }
	// TODO 
	// enforce any specific checks for Map.
    return true;
}


// for map, -c is the column upon which to operate
bool ContextMap::handle_c()
{
    if ((_i+1) < _argc) {
        setColumn(atoi(_argv[_i + 1]));
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}


// for map, -o is the operation to apply to the column (-c)
bool ContextMap::handle_o()
{
    if ((_i+1) < _argc) {
        setColumnOperation(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}


// for map, -null is a NULL vakue assigned
// when no overlaps are detected.
bool ContextMap::handle_null()
{
    if ((_i+1) < _argc) {
        setNullValue(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}
