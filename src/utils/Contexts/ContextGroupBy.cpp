/*
 * ContextGroupBy.cpp
 *
 *  Created on: Mar 26, 2014
 *      Author: nek3d
 */


#include "ContextGroupBy.h"

ContextGroupBy::ContextGroupBy()
: _printFullCols(false),
  _ignoreCase(false)
{
	setSortedInput(true);
	_noEnforceCoordSort = true;
	setColumnOpsMethods(true);
	setExplicitBedOutput(true);

	//For columnOps, groupBy has default operation sum but no default column,
	///so we have to clear that.
	//
	_keyListOps->setColumns("");

}

ContextGroupBy::~ContextGroupBy()
{

}


bool ContextGroupBy::parseCmdArgs(int argc, char **argv, int skipFirstArgs)
{
	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		else if ((strcmp(_argv[_i], "-g") == 0) || (strcmp(_argv[_i], "-grp") == 0)) {
			if (!handle_g()) return false;
		}
		else if (strcmp(_argv[_i], "-inheader") == 0) {
			if (!handle_inheader()) return false;
		}
		else if (strcmp(_argv[_i], "-outheader") == 0) {
			if (!handle_outheader()) return false;
		}
		else if (strcmp(_argv[_i], "-header") == 0) {
			if (!handle_header()) return false;
		}
		else if (strcmp(_argv[_i], "-full") == 0) {
			if (!handle_full()) return false;
		}
		else if (strcmp(_argv[_i], "-ignorecase") == 0) {
			if (!handle_ignorecase()) return false;
		}
	}
	return ContextBase::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextGroupBy::isValidState()
{
	// The user was required to have entered one or more columns
	if (_keyListOps->getColumns().empty()) {
		_errorMsg = "***** ERROR: -opCols parameter requires a value.";
		return false;
	}

	//default to stdin
	if (getNumInputFiles() == 0) {
		addInputFile("-");
	}
	//default grouping is cols 1,2,3
	if (_groupStr.empty()) _groupStr = "1,2,3";

	return ContextBase::isValidState();
}


bool ContextGroupBy::handle_g()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -g option given, but columns to group not specified. *****";
		return false;
	}
	_groupStr = _argv[_i+1];
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_inheader()
{
	_inheader = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_outheader() {
	setPrintHeader(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_header() {
	_inheader = true;
	setPrintHeader(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_full() {
	_printFullCols = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextGroupBy::handle_ignorecase() {
	_ignoreCase = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

const string &ContextGroupBy::getDefaultHeader() {
	//groupBy does not support multiple databases.
	FileRecordMgr *frm = _files[0];
	int numFields = frm->getNumFields();
	_defaultHeader.clear();
	ostringstream s;
	for (int i=1; i <= numFields; i++) {
		s << "col_";
		s << i;
		s << "\t";
	}
	_defaultHeader.append(s.str());
	//change last tab into newline
	_defaultHeader[_defaultHeader.size()-1] = '\n';
	return _defaultHeader;
}
