/*
 * Context.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: nek3d
 */

#include "Context.h"

Context::Context()
:
  _program(UNSPECIFIED_PROGRAM),
  _useMergedIntervals(false),
  _genomeFile(NULL),
  _outputFileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE),
  _outputTypeDetermined(false),
  _skipFirstArgs(0),
  _showHelp(false),
  _obeySplits(false),
  _uncompressedBam(false),
  _anyHit(false),
  _noHit(false),
  _writeA(false),
  _writeB(false),
  _leftJoin(false),
  _writeCount(false),
  _writeOverlap(false),
  _writeAllOverlap(false),
  _haveFraction(false),
  _overlapFraction(1E-9),
  _reciprocal(false),
  _sameStrand(false),
  _diffStrand(false),
   _sortedInput(false),
  _printHeader(false),
  _printable(true),
   _explicitBedOutput(false),
  _queryFileIdx(-1),
  _databaseFileIdx(-1),
  _bamHeaderAndRefIdx(-1),
  _maxNumDatabaseFields(0),
  _reportCount(false),
  _maxDistance(0),
  _reportNames(false),
  _reportScores(false)
{
	_validScoreOps.insert("sum");
	_validScoreOps.insert("max");
	_validScoreOps.insert("min");
	_validScoreOps.insert("mean");
	_validScoreOps.insert("mode");
	_validScoreOps.insert("median");
	_validScoreOps.insert("antimode");
	_validScoreOps.insert("collapse");
}

Context::~Context()
{
	if (_genomeFile != NULL) {
		delete _genomeFile;
		_genomeFile = NULL;
	}
}

bool Context::determineOutputType() {
	if (_outputTypeDetermined) {
		return true;
	}
	//test whether output should be BED or BAM.
	//If the user explicitly requested BED, then it's BED.
	//Otherwise, if there are any BAM files in the input,
	//then the output should be BAM.
	if (getExplicitBedOutput() || getQueryFileType() != FileRecordTypeChecker::BAM_FILE_TYPE) {
		setOutputFileType(FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE);
	} else {
		setOutputFileType(FileRecordTypeChecker::BAM_FILE_TYPE);
	}
	_outputTypeDetermined = true;
	return true;

}

void Context::openGenomeFile(const QuickString &genomeFilename)
{
	_genomeFile = new NewGenomeFile(genomeFilename.c_str());
}

void Context::openGenomeFile(const BamTools::RefVector &refVector)
{
	_genomeFile = new NewGenomeFile(refVector);
}

void Context::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
	_argc = argc;
	_argv = argv;
	_skipFirstArgs = skipFirstArgs;

	_argsProcessed.resize(argc - skipFirstArgs, false);

	for (int i=skipFirstArgs; i < argc; i++) {
		if (isUsed(i - skipFirstArgs)) {
			continue;
		}

		if (strcmp(argv[i], "-i") == 0) {
			addInputFile(argv[i+1]);
			markUsed(i - skipFirstArgs);
			i++;
			markUsed(i - skipFirstArgs);
		} else if (strcmp(argv[i], "-g") == 0) {
			openGenomeFile(argv[i+1]);
			markUsed(i - skipFirstArgs);
			i++;
			markUsed(i - skipFirstArgs);
		} else if (strcmp(argv[i], "-h") == 0) {
			setShowHelp(true);
			markUsed(i - skipFirstArgs);
		} else if (strcmp(argv[i], "--help") == 0) {
			setShowHelp(true);
			markUsed(i - skipFirstArgs);
		}
        else if (strcmp(argv[i], "-split") == 0) {
            setObeySplits(true);
            markUsed(i - skipFirstArgs);
        }
		if (strcmp(argv[i], "-a") == 0) {
			addInputFile(argv[i+1]);
			_queryFileIdx = getNumInputFiles() -1;
			markUsed(i - skipFirstArgs);
			i++;
			markUsed(i - skipFirstArgs);
		}
        else if(strcmp(argv[i], "-abam") == 0) {
			addInputFile(argv[i+1]);
			_queryFileIdx = getNumInputFiles() -1;
			markUsed(i - skipFirstArgs);
			i++;
			markUsed(i - skipFirstArgs);
			setInputFileType(_queryFileIdx, FileRecordTypeChecker::BAM_FILE_TYPE);
        }
        else if (strcmp(argv[i], "-b") == 0) {
			addInputFile(argv[i+1]);
			_databaseFileIdx = getNumInputFiles() -1;
			markUsed(i - skipFirstArgs);
			i++;
			markUsed(i - skipFirstArgs);
		} else if (strcmp(argv[i], "-u") == 0) {
            setAnyHit(true);
            markUsed(i - skipFirstArgs);
        } else if(strcmp(argv[i], "-f") == 0) {
            if ((i+1) < argc) {
                setHaveFraction(true);
                setOverlapFraction(atof(argv[i + 1]));
                markUsed(i - skipFirstArgs);
                i++;
                markUsed(i - skipFirstArgs);
            }
        }
        else if(strcmp(argv[i], "-bed") == 0) {
            setExplicitBedOutput(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-wa") == 0) {
            setWriteA(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-wb") == 0) {
            setWriteB(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-wo") == 0) {
            setWriteOverlap(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-wao") == 0) {
            setWriteAllOverlap(true);
            setWriteOverlap(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-c") == 0) {
            setWriteCount(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-r") == 0) {
            setReciprocal(true);
            markUsed(i - skipFirstArgs);
        }
        else if (strcmp(argv[i], "-v") == 0) {
            setNoHit(true);
            markUsed(i - skipFirstArgs);
        }
        else if (strcmp(argv[i], "-s") == 0) {
            setSameStrand(true);
            markUsed(i - skipFirstArgs);
        }
        else if (strcmp(argv[i], "-S") == 0) {
            setDiffStrand(true);
            markUsed(i - skipFirstArgs);
        }
        else if (strcmp(argv[i], "-loj") == 0) {
            setLeftJoin(true);
            markUsed(i - skipFirstArgs);
        }
        else if (strcmp(argv[i], "-ubam") == 0) {
            setUncompressedBam(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-sorted") == 0) {
            setSortedInput(true);
            markUsed(i - skipFirstArgs);
        }
        else if(strcmp(argv[i], "-header") == 0) {
            setPrintHeader(true);
            markUsed(i - skipFirstArgs);
        }
	}
}

bool Context::isValidState()
{
	if (!Context::cmdArgsValid()) {
		return false;
	}
	if (getAnyHit() && getNoHit()) {
		_errorMsg = "Error: request either -u for anyHit OR -v for noHit, not both.";
		return false;
	}
	if (getWriteCount()) {
		if (getAnyHit()) {
			_errorMsg = "Error: request either -c for writeCount OR -u for anyHit, not both.";
			return false;
		}  else if (getWriteB()) {
			_errorMsg = "Error: request either -c for writeCount OR -wb for writeB, not both.";
			return false;
		} else if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			_errorMsg = "Error: writeCount option is not valid with BAM query input, unless bed output is specified with -bed option.";
			return false;
		}
	}
	if (getWriteOverlap()) {

		if (getWriteA()) {
			_errorMsg = "Error: request either -wa for writeA OR -wo for writeOverlap, not both.";
			return false;
		} else if (getWriteB()) {
			_errorMsg = "Error: request either -wb for writeB OR -wo for writeOverlap, not both.";
			return false;
		}  else if (getWriteCount()) {
			_errorMsg = "Error: request either -c for writeCount OR -wo for writeOverlap, not both.";
			return false;
		} else if (getAnyHit()) {
			_errorMsg = "Error: request either -u for anyHit OR -wo for writeOverlap, not both.";
			return false;
		} else if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			_errorMsg = "Error: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option.";
			return false;
		}
	}
	if (getWriteB() || getLeftJoin()) {
		if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			cerr << endl << "*****" << endl << "*****WARNING: -wb and -loj are ignored with bam input, unless bed output is specified with -bed option." << endl << "*****" << endl;
		}
	}
	if (getSameStrand() && getDiffStrand()) {
		_errorMsg = "Error: request -s for sameStrand, or -S for diffStrand, not both.";
		return false;
	}

	if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && getPrintHeader()) {
		cerr << endl << "*****" << endl << "*****WARNING: -header option is not valid for BAM input." << endl << "*****" << endl;
		setPrintHeader(false);
	}
	if (getAnyHit() || getNoHit() || getWriteCount()) {
		setPrintable(false);
	}
	return _inputFiles.size() > 0;
}

bool Context::cmdArgsValid()
{
	bool retval = true;
	for (int i = _skipFirstArgs; i < _argc; i++) {
		if (!isUsed(i - _skipFirstArgs)) {
			_errorMsg += "\nERROR. Unrecognized argument: ";
			_errorMsg += _argv[i];
			retval = false;
		}
	}
	return retval;
}

int Context::getBamHeaderAndRefIdx() {
	if (_bamHeaderAndRefIdx != -1) {
		//already found which BAM file to usefor the header
		return _bamHeaderAndRefIdx;
	}
	if (_inputFiles[_queryFileIdx]._fileType == FileRecordTypeChecker::BAM_FILE_TYPE) {
		_bamHeaderAndRefIdx = _queryFileIdx;
	} else {
		_bamHeaderAndRefIdx = _databaseFileIdx;
	}
	return _bamHeaderAndRefIdx;
}




