/*
 * ContextFisher.cpp
 *
 */

#include "ContextFisher.h"

ContextFisher::ContextFisher() {
	setSortedInput(true);
	setUseMergedIntervals(false);
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
//        if (strcmp(_argv[_i], "-g") == 0) {
//              if (!handle_g()) return false;
//        }
//        if(strcmp(_argv[_i], "-m") == 0) {
//            markUsed(_i - _skipFirstArgs);
//	        setUseMergedIntervals(true);
//        }
	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextFisher::isValidState()
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
    if (_genomeFile == NULL){
        _errorMsg = "\nERROR*****: specify -g genome file*****\n";
        return false;
    }
	//column operations not allowed with BAM input
	if (hasColumnOpsMethods() &&
			getFile(0)->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		_errorMsg = "\n***** ERROR: stranded merge not supported for VCF files. *****";
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

