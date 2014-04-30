/*
 * ContextJaccard.cpp
 *
 *  Created on: Apr 24, 2014
 *      Author: nek3d
 */

#include "ContextJaccard.h"

ContextJaccard::ContextJaccard() {
	setSortedInput(true);
	setUseMergedIntervals(true);

}

ContextJaccard::~ContextJaccard() {

}

bool ContextJaccard::parseCmdArgs(int argc, char **argv, int skipFirstArgs)
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
		else if (strcmp(_argv[_i], "-s") == 0) {
			if (!handle_s()) return false;
		}
		else if (strcmp(_argv[_i], "-S") == 0) {
			if (!handle_S()) return false;
		}
	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextJaccard::isValidState()
{
	if (!ContextIntersect::isValidState()) {
		return false;
	}
	// Tests for stranded merge
	//
	if (_desiredStrand != FileRecordMergeMgr::ANY_STRAND) { // requested stranded merge
		for (int i=0; i < getNumInputFiles(); i++) {
			// make sure file has strand.
			if (!getFile(i)->recordsHaveStrand()) {
				_errorMsg = "\n***** ERROR: stranded merge requested, but input file ";
				_errorMsg  += getInputFileName(i);
				_errorMsg  += " does not have strands. *****";
				return false;
			}
			//make sure file is not VCF.
			if (getFile(1)->getFileType() == FileRecordTypeChecker::VCF_FILE_TYPE) {
				_errorMsg = "\n***** ERROR: stranded merge not supported for VCF file ";
				_errorMsg += getInputFileName(i);
				_errorMsg += ". *****";
				return false;
			}
		}
	}
	//column operations not allowed with BAM input
	if (hasColumnOpsMethods() &&
			getFile(0)->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		_errorMsg = "\n***** ERROR: stranded merge not supported for VCF files. *****";
		return false;
	}
	return true;
}

bool ContextJaccard::handle_s() {
	setSameStrand(true);
	_desiredStrand = FileRecordMergeMgr::SAME_STRAND_EITHER;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextJaccard::handle_S() {
    if ((_i+1) < _argc) {
    	bool validChar = false;
    	if (_argv[_i+1][0] == '+') {
			_desiredStrand = FileRecordMergeMgr::SAME_STRAND_FORWARD;
			validChar = true;
    	} else if (_argv[_i+1][0] == '-') {
    		validChar = true;
			_desiredStrand = FileRecordMergeMgr::SAME_STRAND_REVERSE;
    	}
    	if (validChar) {
			markUsed(_i - _skipFirstArgs);
			_i++;
			markUsed(_i - _skipFirstArgs);
			setSameStrand(true);
			return true;
    	}
    }
	_errorMsg = "\n***** ERROR: -S option must be followed by + or -. *****";
	return false;
}
