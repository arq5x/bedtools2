/*
 * ContextBase.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: nek3d
 */

#include "ContextBase.h"
#include <unistd.h>
#include <sys/types.h>

ContextBase::ContextBase()
:
  _program(UNSPECIFIED_PROGRAM),
  _allFilesOpened(false),
  _useMergedIntervals(false),
  _genomeFile(NULL),
  _outputFileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE),
  _outputTypeDetermined(false),
  _skipFirstArgs(0),
  _showHelp(false),
  _obeySplits(false),
  _uncompressedBam(false),
  _useBufferedOutput(true),
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
  _useFullBamTags(false),
  _reportCount(false),
  _maxDistance(0),
  _reportNames(false),
  _reportScores(false),
  _numOutputRecords(0),
  _hasConstantSeed(false),
  _seed(0),
  _forwardOnly(false),
  _reverseOnly(false)
{
	_programNames["intersect"] = INTERSECT;
	_programNames["sample"] = SAMPLE;
	_programNames["map"] = MAP;

	_validScoreOps.insert("sum");
	_validScoreOps.insert("max");
	_validScoreOps.insert("min");
	_validScoreOps.insert("mean");
	_validScoreOps.insert("mode");
	_validScoreOps.insert("median");
	_validScoreOps.insert("antimode");
	_validScoreOps.insert("collapse");
}

ContextBase::~ContextBase()
{
	delete _genomeFile;
	_genomeFile = NULL;

	//close all files and delete FRM objects.
	for (int i=0; i < (int)_files.size(); i++) {
		_files[i]->close();
		delete _files[i];
		_files[i] = NULL;
	}
}

bool ContextBase::determineOutputType() {
	if (_outputTypeDetermined) {
		return true;
	}
	//test whether output should be BED or BAM.
	//If the user explicitly requested BED, then it's BED.
	if (getExplicitBedOutput()) {
		setOutputFileType(FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE);
		_outputTypeDetermined = true;
		return true;
	}

	//Otherwise, if there are any BAM files in the input,
	//then the output should be BAM.
	for (_i = 0; _i < (int)_files.size(); _i++) {
		if (_files[_i]->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
			setOutputFileType(FileRecordTypeChecker::BAM_FILE_TYPE);
			_bamHeaderAndRefIdx = _i;
			_outputTypeDetermined = true;
			return true;
		}
	}

	//Okay, it's bed.
	setOutputFileType(FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE);
	_outputTypeDetermined = true;
	return true;


}

void ContextBase::openGenomeFile(const QuickString &genomeFilename)
{
	_genomeFile = new NewGenomeFile(genomeFilename.c_str());
}

void ContextBase::openGenomeFile(const BamTools::RefVector &refVector)
{
	_genomeFile = new NewGenomeFile(refVector);
}

bool ContextBase::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
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

		if (strcmp(_argv[_i], "-i") == 0) {
			if (!handle_i()) return false;
		}
		else if (strcmp(_argv[_i], "-g") == 0) {
			if (!handle_g()) return false;
		}
		else if ((strcmp(_argv[_i], "-h") == 0) || (strcmp(_argv[_i], "--help") == 0)) {
			if (!handle_h()) return false;
		}
		else if (strcmp(_argv[_i], "-split") == 0) {
			if (!handle_split()) return false;
		}
        else if (strcmp(_argv[_i], "-bed") == 0) {
			if (!handle_bed()) return false;
       }
        else if (strcmp(_argv[_i], "-ubam") == 0) {
			if (!handle_ubam()) return false;
        }
        else if (strcmp(_argv[_i], "-fbam") == 0) {
			if (!handle_fbam()) return false;
        }
        else if(strcmp(_argv[_i], "-sorted") == 0) {
			if (!handle_sorted()) return false;
        }
        else if (strcmp(_argv[_i], "-nobuf") == 0) {
			if (!handle_nobuf()) return false;
        }
        else if (strcmp(_argv[_i], "-header") == 0) {
			if (!handle_header()) return false;
        }
        else if (strcmp(_argv[_i], "-n") == 0) {
			if (!handle_n()) return false;
        }
        else if (strcmp(_argv[_i], "-seed") == 0) {
			if (!handle_seed()) return false;
        }
	}
	return true;
}

bool ContextBase::isValidState()
{
	if (!openFiles()) {
		return false;
	}
	if (!cmdArgsValid()) {
		return false;
	}
	if (!determineOutputType()) {
		return false;
	}
	return true;
}


bool ContextBase::cmdArgsValid()
{
	bool retval = true;
	for (_i = _skipFirstArgs; _i < _argc; _i++) {
		if (!isUsed(_i - _skipFirstArgs)) {
			_errorMsg += "\n***** ERROR: Unrecognized parameter: ";
			_errorMsg += _argv[_i];
			_errorMsg += " *****";
			retval = false;
		}
	}
	return retval;
}

bool ContextBase::openFiles() {

	//Make a vector of FileRecordMgr objects by going through the vector
	//of filenames and opening each one.
	if (_allFilesOpened) {
		return true;
	}
	_files.resize(_fileNames.size());

	for (int i = 0; i < (int)_fileNames.size(); i++) {
		FileRecordMgr *frm = new FileRecordMgr(_fileNames[i], _sortedInput);
		if (hasGenomeFile()) {
			frm->setGenomeFile(_genomeFile);
		}
		frm->setFullBamFlags(_useFullBamTags);
		if (!frm->open()) {
			return false;
		}
		_files[i] = frm;
	}
	_allFilesOpened = true;
	return true;
}

int ContextBase::getBamHeaderAndRefIdx() {
	if (_bamHeaderAndRefIdx != -1) {
		//already found which BAM file to use for the header
		return _bamHeaderAndRefIdx;
	}
	if (_files[_queryFileIdx]->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		_bamHeaderAndRefIdx = _queryFileIdx;
	} else {
		_bamHeaderAndRefIdx = _databaseFileIdx;
	}
	return _bamHeaderAndRefIdx;
}

int ContextBase::getUnspecifiedSeed()
{
	// thanks to Rob Long for the tip.
	_seed = (unsigned)time(0)+(unsigned)getpid();
	srand(_seed);
	return _seed;
}

bool ContextBase::handle_bed()
{
	setExplicitBedOutput(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_fbam()
{
	setUseFullBamTags(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_g()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -g option given, but no genome file specified. *****";
		return false;
	}
	openGenomeFile(_argv[_i+1]);
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_h()
{
	setShowHelp(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_header()
{
	setPrintHeader(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_i()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -i option given, but no input file specified. *****";
		return false;
	}
	addInputFile(_argv[_i+1]);
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_n()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -n option given, but no number of output records specified. *****";
		return false;
	}
	setNumOutputRecords(atoi(_argv[_i + 1]));
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_nobuf()
{
	setUseBufferedOutput(false);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_seed()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -seed option given, but no seed specified. *****";
		return false;
	}
	_hasConstantSeed = true;
	_seed  = atoi(_argv[_i+1]);
	srand(_seed);
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_split()
{
    setObeySplits(true);
    markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_sorted()
{
	setSortedInput(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_ubam()
{
    setUncompressedBam(true);
    markUsed(_i - _skipFirstArgs);
	return true;
}
