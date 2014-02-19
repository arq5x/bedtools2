/*
 * ContextMap.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: arq5x
 */

#include "ContextMap.h"

ContextMap::ContextMap()
: _delimStr(",")
{
	// map requires sorted input
	setSortedInput(true);
	setLeftJoin(true);

	// default to BED score column
	setColumns("5");
	// default to "sum"
	setOperations("sum");
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
        else if (strcmp(_argv[_i], "-delim") == 0) {
			if (!handle_delim()) return false;
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
        cerr << endl << "*****" << endl
             << "***** ERROR: BAM database file not currently supported for the map tool." 
             << endl;
        exit(1);
    }


	//get the strings from context containing the comma-delimited lists of columns
	//and operations. Split both of these into vectors. Get the operation code
	//for each operation string. Finally, make a vector of pairs, where the first
	//member of each pair is a column number, and the second member is the code for the
	//operation to perform on that column.

	vector<QuickString> columnsVec;
	vector<QuickString> opsVec;
	int numCols = Tokenize(_columns, columnsVec, ',');
	int numOps = Tokenize(_operations, opsVec, ',');

	if (numOps < 1 || numCols < 1) {
		 cerr << endl << "*****" << endl
		             << "***** ERROR: There must be at least one column and at least one operation named." << endl;
		 return false;
	}
	if (numOps > 1 && numCols != numOps) {
		 cerr << endl << "*****" << endl
		             << "***** ERROR: There are " << numCols <<" columns given, but there are " << numOps << " operations. " << endl;
		cerr << "\tPlease provide either a single operation that will be applied to all listed columns, " << endl;
		cerr << "\tor an operation for each column." << endl;
		return false;
	}
	KeyListOps keyListOps;
	for (int i=0; i < (int)columnsVec.size(); i++) {
		int col = str2chrPos(columnsVec[i]);

		//check that the column number is valid
		if (col < 1 || col > getDatabaseFile()->getNumFields()) {
			 cerr << endl << "*****" << endl  << "***** ERROR: Requested column " << col << ", but database file "
					 << getDatabaseFileName() << " only has fields 1 - " << getDatabaseFile()->getNumFields() << "." << endl;
			 return false;
		}
		const QuickString &operation = opsVec.size() > 1 ? opsVec[i] : opsVec[0];
		KeyListOps::OP_TYPES opCode = keyListOps.getOpCode(operation);
		if (opCode == KeyListOps::INVALID) {
			cerr << endl << "*****" << endl
								 << "***** ERROR: " << operation << " is not a valid operation. " << endl;
			return false;
		}
		_colOps.push_back(pair<int, KeyListOps::OP_TYPES>(col, opCode));
	}
    return true;
}


// for map, -c is the string of columns upon which to operate
bool ContextMap::handle_c()
{
    if ((_i+1) < _argc) {
        setColumns(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}


// for map, -o is the string of operations to apply to the columns (-c)
bool ContextMap::handle_o()
{
    if ((_i+1) < _argc) {
        setOperations(_argv[_i + 1]);
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

bool ContextMap::handle_delim()
{
    if ((_i+1) < _argc) {
    	_delimStr = _argv[_i + 1];
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}
