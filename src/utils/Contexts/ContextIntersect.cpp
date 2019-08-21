/*
 * ContextIntersect.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: nek3d
 */

#include "ContextIntersect.h"

ContextIntersect::ContextIntersect()
{

}

ContextIntersect::~ContextIntersect()
{

}


bool ContextIntersect::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		if (strcmp(_argv[_i], "-a") == 0) {
			if (!handle_a()) return false;
		}
        else if (strcmp(_argv[_i], "-abam") == 0) {
			if (!handle_abam()) return false;
        }
        else if (strcmp(_argv[_i], "-b") == 0) {
			if (!handle_b()) return false;
		}
        else if (strcmp(_argv[_i], "-names") == 0) {
			if (!handle_names()) return false;
		}
        else if (strcmp(_argv[_i], "-filenames") == 0) {
			if (!handle_filenames()) return false;
		}
        else if (strcmp(_argv[_i], "-u") == 0) {
			if (!handle_u()) return false;
        }
        else if (strcmp(_argv[_i], "-f") == 0) {
			if (!handle_f()) return false;
        }
        else if (strcmp(_argv[_i], "-F") == 0) {
			if (!handle_F()) return false;
        }
        else if (strcmp(_argv[_i], "-wa") == 0) {
			if (!handle_wa()) return false;
        }
        else if (strcmp(_argv[_i], "-wao") == 0) {
			if (!handle_wao()) return false;
        }
        else if (strcmp(_argv[_i], "-wb") == 0) {
			if (!handle_wb()) return false;
        }
        else if (strcmp(_argv[_i], "-wo") == 0) {
			if (!handle_wo()) return false;
        }
        else if (strcmp(_argv[_i], "-c") == 0) {
			if (!handle_c()) return false;
        }
        else if (strcmp(_argv[_i], "-C") == 0) {
			if (!handle_C()) return false;
        }
        else if(strcmp(_argv[_i], "-r") == 0) {
			if (!handle_r()) return false;
        }
        else if(strcmp(_argv[_i], "-e") == 0) {
			if (!handle_e()) return false;
        }
        else if (strcmp(_argv[_i], "-v") == 0) {
			if (!handle_v()) return false;
        }
        else if (strcmp(_argv[_i], "-s") == 0) {
			if (!handle_s()) return false;
        }
        else if (strcmp(_argv[_i], "-S") == 0) {
			if (!handle_S()) return false;
        }
        else if (strcmp(_argv[_i], "-loj") == 0) {
			if (!handle_loj()) return false;
        }
	}
	return ContextBase::parseCmdArgs(argc, argv, _skipFirstArgs);
}



bool ContextIntersect::isValidState()
{
	if (!ContextBase::isValidState()) {
		return false;
	}

	if (_queryFileIdx == -1 || _dbFileIdxs.size() == -0) {
		_errorMsg = "\n***** ERROR: query and database files not specified. *****";
		return false;
	}

	if (getAnyHit() && getNoHit()) {
		_errorMsg = "\n***** ERROR: request either -u for anyHit OR -v for noHit, not both. *****";
		return false;
	}
	if (getWriteCount()) {
		if (getAnyHit()) {
			_errorMsg = "\n***** ERROR: request either -c for writeCount OR -u for anyHit, not both. *****";
			return false;
		}  else if (getWriteB()) {
			_errorMsg = "\n***** ERROR: request either -c for writeCount OR -wb for writeB, not both. *****";
			return false;
		} else if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			_errorMsg = "\n***** ERROR: writeCount option is not valid with BAM query input, unless bed output is specified with -bed option. *****";
			return false;
		}
	}

	if (_haveFractionA && (_overlapFractionA <= 0.0 || _overlapFractionA > 1.0)) {
		_errorMsg = "\n***** ERROR: -f must be in the range (0.0, 1.0]. *****";
		return false;
	}
	if (_haveFractionB && (_overlapFractionB <= 0.0 || _overlapFractionB > 1.0)) {
		_errorMsg = "\n***** ERROR: -F must be in the range (0.0, 1.0]. *****";
		return false;
	}
	if (_haveFractionB && !_haveFractionA && _eitherFraction) {
		_errorMsg = "\n***** ERROR: -e must be used with -f. *****";
		return false;
	}
	if (_haveFractionB && _reciprocalFraction) {
		_errorMsg = "\n***** ERROR: -r must be used solely with -f. *****";
		return false;
	}
	if (_haveFractionA && _haveFractionB && _reciprocalFraction) {
		_errorMsg = "\n***** ERROR: -r must be used solely with -f. *****";
		return false;
	}
	if (_writeCount && (_reportDBfileNames || _reportDBnameTags)) {
		_errorMsg = "\n***** ERROR: -c cannot be used with -filenames or -names. Use -C. *****";
		return false;
	}
	// if -f and -r are given, then we need to set _overlapFractionB to be the same as _overlapFractionA
	if (_haveFractionA && _reciprocalFraction) {
		setOverlapFractionB(_overlapFractionA);
	}

	if (getUseDBnameTags() && _dbNameTags.size() != _dbFileIdxs.size()) {
		_errorMsg = "\n***** ERROR: Number of database name tags given does not match number of databases. *****";
		return false;

	}


	if (getWriteOverlap()) {

		if (getWriteA()) {
			_errorMsg = "\n***** ERROR: request either -wa for writeA OR -wo for writeOverlap, not both. *****";
			return false;
		} else if (getWriteB()) {
			_errorMsg = "\n***** ERROR: request either -wb for writeB OR -wo for writeOverlap, not both. *****";
			return false;
		}  else if (getWriteCount()) {
			_errorMsg = "\n***** ERROR: request either -c for writeCount OR -wo for writeOverlap, not both. *****";
			return false;
		} else if (getAnyHit()) {
			_errorMsg = "\n***** ERROR: request either -u for anyHit OR -wo for writeOverlap, not both. *****";
			return false;
		} else if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			_errorMsg = "\n***** ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option. *****";
			return false;
		}
	}
	if (getWriteB() || getLeftJoin()) {
		if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			cerr << endl << "*****" << endl << "*****WARNING: -wb and -loj are ignored with bam input, unless bed output is specified with -bed option." << endl << "*****" << endl;
		}
	}
	if (getSameStrand() && getDiffStrand()) {
		_errorMsg = "\n***** ERROR: request -s for sameStrand, or -S for diffStrand, not both. *****";
		return false;
	}
	if ((getSameStrand() || getDiffStrand()) && !strandedToolSupported()) {
		//error msg set within strandedToolSupported method.
		return false;
	}

	if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && getPrintHeader()) {
		cerr << endl << "*****" << endl << "*****WARNING: -header option is not valid for BAM input." << endl << "*****" << endl;
		setPrintHeader(false);
	}
	if (getAnyHit() || getNoHit() || getWriteCount()) {
		setPrintable(false);
	}

	if (getNoHit() || getWriteCount() || getWriteOverlap() 
		|| getWriteAllOverlap() || getLeftJoin())
	{
		setRunToQueryEnd(true);
	}

	if (_files.size()  < 2 ) {
		return false;
	}
	return true;
}

bool ContextIntersect::determineOutputType() {
	if (_outputTypeDetermined) {
		return true;
	}

	//determine the maximum number of database fields.
	for (int i=0; i < getNumDatabaseFiles(); i++) {
		int numFields = getDatabaseFile(i)->getNumFields();
		if ( numFields > _maxNumDatabaseFields) {
			_maxNumDatabaseFields = numFields;
		}
	}
	return ContextBase::determineOutputType();
}

bool ContextIntersect::handle_a()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -a option given, but no query file specified. *****";
		return false;
	}

	addInputFile(_argv[_i+1]);
	_queryFileIdx = getNumInputFiles() -1;
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextIntersect::handle_abam()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -abam option given, but no query BAM file specified. *****";
		return false;
	}
	addInputFile(_argv[_i+1]);
	_queryFileIdx = getNumInputFiles() -1;
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextIntersect::handle_b()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -b option given, but no database file specified. *****";
		return false;
	}

	do {
		addInputFile(_argv[_i+1]);
		int fileId = getNumInputFiles() -1;
		_dbFileIdxs.push_back(fileId);
		_fileIdsToDbIdxs[fileId] = _dbFileIdxs.size() -1;
		markUsed(_i - _skipFirstArgs);
		_i++;
		markUsed(_i - _skipFirstArgs);
	} while (_argc > _i+1 && _argv[_i+1][0] != '-');
	return true;
}


bool ContextIntersect::handle_names()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -b option given, but no database names specified. *****";
		return false;
	}

	do {
		addDatabaseNameTag(_argv[_i+1]);
		markUsed(_i - _skipFirstArgs);
		_i++;
		markUsed(_i - _skipFirstArgs);
	} while (_argc > _i+1 && _argv[_i+1][0] != '-');
	setUseDBnameTags(true);
	return true;
}

bool ContextIntersect::handle_filenames()
{
	markUsed(_i - _skipFirstArgs);
	setUseDBfileNames(true);
	return true;
}

bool ContextIntersect::handle_c()
{
    setWriteCount(true);
    markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_C()
{
    setWriteCountsPerDatabase(true);
    markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_f()
{
    if ((_i+1) < _argc) {
        setHaveFractionA(true);
        setOverlapFractionA(atof(_argv[_i + 1]));
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}

bool ContextIntersect::handle_F()
{
    if ((_i+1) < _argc) {
        setHaveFractionB(true);
        setOverlapFractionB(atof(_argv[_i + 1]));
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}

bool ContextIntersect::handle_loj()
{
    setLeftJoin(true);
    markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_r()
{
	setReciprocalFraction(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_e()
{
	setEitherFraction(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_s()
{
	setSameStrand(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_S()
{
	setDiffStrand(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_u()
{
    setAnyHit(true);
    markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_v()
{
	setNoHit(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_wa()
{
	setWriteA(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_wao()
{
	setWriteAllOverlap(true);
	setWriteOverlap(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_wb()
{
    setWriteB(true);
     markUsed(_i - _skipFirstArgs);
     return true;
}

bool ContextIntersect::handle_wo()
{
	setWriteOverlap(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}
