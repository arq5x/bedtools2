/*
 * ContextBase.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: nek3d
 */

#include "ContextBase.h"
#include "Random.h"
#include <unistd.h>
#include <sys/types.h>
#include <cctype>

ContextBase::ContextBase()
:
  _program(UNSPECIFIED_PROGRAM),
  _allFilesOpened(false),
  _genomeFile(NULL),
  _outputFileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE),
  _outputTypeDetermined(false),
  _skipFirstArgs(0),
  _showHelp(false),
  _obeySplits(false),
  _uncompressedBam(false),
  _useBufferedOutput(true),
  _ioBufSize(0),
  _anyHit(false),
  _noHit(false),
  _writeA(false),
  _writeB(false),
  _leftJoin(false),
  _writeCount(false),
  _writeCountsPerDatabase(false),
  _writeOverlap(false),
  _writeAllOverlap(false),
  _haveFractionA(false),
  _haveFractionB(false),
  _overlapFractionA(0.0),
  _overlapFractionB(0.0),
  _reciprocalFraction(false),
  _eitherFraction(false),
  _sameStrand(false),
  _diffStrand(false),
  _sortedInput(false),
  _sortOutput(false),
  _reportDBnameTags(false),
  _reportDBfileNames(false),
  _printHeader(false),
  _printable(true),
  _explicitBedOutput(false),
  _queryFileIdx(-1),
  _bamHeaderAndRefIdx(-1),
  _maxNumDatabaseFields(0),
  _useFullBamTags(false),
  _numOutputRecords(0),
  _hasConstantSeed(false),
  _seed(0),
  _forwardOnly(false),
  _reverseOnly(false),
  _nameCheckDisabled(false),
  _hasColumnOpsMethods(false),
  _keyListOps(NULL),
  _desiredStrand(FileRecordMergeMgr::ANY_STRAND),
  _maxDistance(0),
  _useMergedIntervals(false),
  _reportPrecision(-1),
  _splitBlockInfo(NULL),
  _allFilesHaveChrInChromNames(UNTESTED),
  _allFileHaveLeadingZeroInChromNames(UNTESTED),
  _noEnforceCoordSort(false),
  _inheader(false),
  _nameConventionWarningTripped(false)

{
	_programNames["intersect"] = INTERSECT;
	_programNames["sample"] = SAMPLE;
	_programNames["map"] = MAP;
	_programNames["merge"] = MERGE;
	_programNames["closest"] = CLOSEST;
	_programNames["subtract"] = SUBTRACT;
	_programNames["jaccard"] = JACCARD;
	_programNames["spacing"] = SPACING;
	_programNames["fisher"] = FISHER;
	_programNames["sample"] = SAMPLE;
	_programNames["coverage"] = COVERAGE;
	_programNames["complement"] = COMPLEMENT;
	_programNames["groupby"] = GROUP_BY;



	if (hasColumnOpsMethods()) {
		_keyListOps = new KeyListOps();
	}
}

ContextBase::~ContextBase()
{
	delete _genomeFile;
	_genomeFile = NULL;
	delete _splitBlockInfo;
	_splitBlockInfo = NULL;

	if (_nameConventionWarningTripped) {
		cerr << _nameConventionWarningMsg << endl;
	}

	//close all files and delete FRM objects.
	for (int i=0; i < (int)_files.size(); i++) {
		_files[i]->close();
		delete _files[i];
		_files[i] = NULL;
	}
	if (hasColumnOpsMethods()) {
		delete _keyListOps;
		_keyListOps = NULL;
	}
}

bool ContextBase::errorEncountered() {
	// just a subcommand was given with no options.
	if (_argc == 1) 
	{
		return true;
	} 
	return !_errorMsg.empty() || getShowHelp();
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

	//Otherwise, if the input is BAM, then the output is BAM
	int fileIdx = hasIntersectMethods() ? _queryFileIdx : 0;
	if (_files[fileIdx]->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		setOutputFileType(FileRecordTypeChecker::BAM_FILE_TYPE);
		_isCram = _files[fileIdx]->isCram();
		return true;
	}

	//Okay, it's bed.
	setOutputFileType(FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE);
	_outputTypeDetermined = true;
	return true;


}

void ContextBase::openGenomeFile(const string &genomeFilename)
{
	_genomeFile = new NewGenomeFile(genomeFilename.c_str());
}

void ContextBase::openGenomeFile(const BamTools::RefVector &refVector)
{
	_genomeFile = new NewGenomeFile(refVector);
}

bool ContextBase::testCmdArgs(int argc, char **argv) {
	_argc = argc;
	_argv = argv;
	_skipFirstArgs = 1;
	_origProgramName = argv[0];
	setProgram(_programNames[_origProgramName]);
	_argsProcessed.resize(_argc - _skipFirstArgs, false);

	if (!parseCmdArgs(argc, argv, 1) || getShowHelp() || !isValidState()) {
		return false;
	}
	return true;
}

bool ContextBase::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {


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
        else if (strcmp(_argv[_i], "-iobuf") == 0) {
			if (!handle_iobuf()) return false;
        }
        else if (strcmp(_argv[_i], "-prec") == 0) {
			if (!handle_prec()) return false;
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
        else if ((strcmp(_argv[_i], "-o") == 0) || (strcmp(_argv[_i], "-ops") == 0)) {
			if (!handle_o()) return false;
        }
        else if ((strcmp(_argv[_i], "-c") == 0) || (strcmp(_argv[_i], "-opCols") == 0)) {
			if (!handle_c()) return false;
        }
        else if (strcmp(_argv[_i], "-null") == 0) {
			if (!handle_null()) return false;
        }
        else if (strcmp(_argv[_i], "-delim") == 0) {
			if (!handle_delim()) return false;
        }
        else if (strcmp(_argv[_i], "-sortout") == 0) {
			if (!handle_sortout()) return false;
        }
        else if (strcmp(_argv[_i], "-nonamecheck") == 0) {
			if (!handle_nonamecheck()) return false;
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
	if (_program != GROUP_BY && 
		_files[0]->getRecordType() == FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE) 
	{
		_errorMsg = "ERROR: file ";
		_errorMsg.append(_files[0]->getFileName());
		_errorMsg.append(" has non positional records, which are only valid for \n");
		_errorMsg.append(" the groupBy tool. Perhaps you are using a header");
		_errorMsg.append(" line(s) that starts with \n");
		_errorMsg.append(" something other than \"#\", \"chrom\", or \"chr\" (any case)?");
		return false;
	}
	if (getObeySplits()) {
		_splitBlockInfo = new BlockMgr(_overlapFractionA, _overlapFractionB, _reciprocalFraction);
	}
	if (hasColumnOpsMethods()) {

		if (hasIntersectMethods()) {
			for (int i=0; i < (int)_dbFileIdxs.size(); i++) {
				FileRecordMgr *dbFile = getFile(_dbFileIdxs[i]);
				_keyListOps->setDBfileType(dbFile->getFileType());
				if (!_keyListOps->isValidColumnOps(dbFile)) {
					return false;
				}
			}
		} else {
			FileRecordMgr *dbFile = getFile(0);
			_keyListOps->setDBfileType(dbFile->getFileType());
			if (!_keyListOps->isValidColumnOps(dbFile)) {
				return false;
			}
		}
		//if user specified a precision, pass it to
		//keyList ops
		if (_reportPrecision != -1) {
			_keyListOps->setPrecision(_reportPrecision);
		}
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

	if (_fileNames.size() == 0) {
		//No input was specified. Error and exit.
		_errorMsg += "\n***** ERROR: No input file given. Exiting. *****";
		return false;
	}

	_files.resize(_fileNames.size());

	for (int i = 0; i < (int)_fileNames.size(); i++) {
		FileRecordMgr *frm = getNewFRM(_fileNames[i], i);
		if (hasGenomeFile()) {
			frm->setGenomeFile(_genomeFile);
		}
		//If we're going to do column operations, and an input file
		// is BAM, we'll need the full flags.
		if (hasColumnOpsMethods()) {
			setUseFullBamTags(true);
		}
		frm->setFullBamFlags(_useFullBamTags);
		frm->setIsSorted(_sortedInput);
		frm->setIoBufSize(_ioBufSize);
		frm->setIsGroupBy(_program == GROUP_BY);
		if (!frm->open(_inheader)) {
			return false;
		}
		if (_noEnforceCoordSort) {
			frm->setNoEnforceCoordSort(true);
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
	if (hasIntersectMethods()) {
		if (_files[_queryFileIdx]->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
			_bamHeaderAndRefIdx = _queryFileIdx;
		} else {
			_bamHeaderAndRefIdx = _dbFileIdxs[0];
		}
		return _bamHeaderAndRefIdx;
	}
	if (_files[0]->getFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		_bamHeaderAndRefIdx = 0;
		return _bamHeaderAndRefIdx;
	}

	return _bamHeaderAndRefIdx;
}

int ContextBase::getUnspecifiedSeed()
{
	// thanks to Rob Long for the tip.
	_seed = (unsigned)time(0)+(unsigned)getpid();
	rand_set_seed(_seed);
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

bool ContextBase::handle_iobuf()
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -iobuf option given, but size of input buffer not specified. *****";
		return false;
	}
	if (!parseIoBufSize(_argv[_i + 1])) return false;
	markUsed(_i - _skipFirstArgs);
	_i++;
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
	rand_set_seed(_seed);
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


// Methods specific to column operations.
// for col ops, -c is the string of columns upon which to operate
bool ContextBase::handle_c()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
        _keyListOps->setColumns(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
    return false;
}


// for col ops, -o is the string of operations to apply to the columns (-c)
bool ContextBase::handle_o()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
    	 _keyListOps->setOperations(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
    return false;
}

bool ContextBase::handle_prec()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
    	int prec = atoi(_argv[_i + 1]);
    	if (prec < 1) {
    		_errorMsg += "\n***** ERROR: -prec must be followed by a positive integer. Exiting. *****";
    		return false;
    	}
    	 _reportPrecision = prec;
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
	_errorMsg += "\n***** ERROR: -prec must be followed by a positive integer. Exiting. *****";
    return false;
}



// for col ops, -null is a NULL value assigned
// when no overlaps are detected.
bool ContextBase::handle_null()
{
	if (!hasColumnOpsMethods()) {
		return false;
	}
    if ((_i+1) < _argc) {
    	 _keyListOps->setNullValue(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
        return true;
    }
    return false;
}

//for col ops, delimStr will appear between each item in
//a collapsed but delimited list.
bool ContextBase::handle_delim()
{
	if (!hasColumnOpsMethods()) {
		_errorMsg = "\n***** ERROR: Can't set delimiter for tools without column operations. Exiting. *****";
		return false;
	}
    if ((_i+1) < _argc) {
    	 _keyListOps->setDelimStr(_argv[_i + 1]);
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}

bool ContextBase::handle_sortout()
{
	setSortOutput(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextBase::handle_nonamecheck()
{
	setNameCheckDisabled(true);
	markUsed(_i - _skipFirstArgs);
	return true;
}

void ContextBase::setColumnOpsMethods(bool val)
{
	if (val && !_hasColumnOpsMethods) {
		//was off, but we're turning it on.
		_keyListOps = new KeyListOps();
	}
	_hasColumnOpsMethods = val;
}

const string &ContextBase::getColumnOpsVal(RecordKeyVector &keyList) const {
	if (!hasColumnOpsMethods()) {
		return _nullStr;
	}
	return _keyListOps->getOpVals(keyList);
}

FileRecordMgr *ContextBase::getNewFRM(const string &filename, int fileIdx) {

	if (_useMergedIntervals) {
		FileRecordMergeMgr *frm = new FileRecordMergeMgr(filename);
		frm->setStrandType(_desiredStrand);
		frm->setMaxDistance(_maxDistance);
		frm->setFileIdx(fileIdx);
		return frm;
	} else {
		FileRecordMgr *frm = new FileRecordMgr(filename);
		frm->setFileIdx(fileIdx);
		return frm;
	}
}

bool ContextBase::parseIoBufSize(string bufStr)
{
	char lastChar = bufStr[bufStr.size()-1];
	size_t multiplier = 1;
	if (!isdigit(lastChar)) {
		switch (lastChar) {
		case 'K':
			multiplier = 1 << 10;
			break;
		case 'M':
			multiplier = 1 << 20;
			break;
		case 'G':
			multiplier = 1 << 30;
			break;
		default:
			_errorMsg = "\n***** ERROR: Unrecognized memory buffer size suffix \'";
			_errorMsg += lastChar;
			_errorMsg += "\' given. *****";
			return false;
			break;
		}
		//lop off suffix character
		bufStr.resize(bufStr.size()-1);
	}
	if (!isNumeric(bufStr)) {
		_errorMsg = "\n***** ERROR: argument passed to -iobuf is not numeric. *****";
		return false;
	}
	_ioBufSize = (size_t)str2chrPos(bufStr) * multiplier;
	if (_ioBufSize < MIN_ALLOWED_BUF_SIZE) {
		_errorMsg = "\n***** ERROR: specified buffer size is too small. *****";
		return false;
	}
	return true;
}

void ContextBase::testNameConventions(const Record *record) {
	//Do nothing if using the -nonamecheck option,
	//warning already given, or record is unmapped BAM record
	if (getNameCheckDisabled() || _nameConventionWarningTripped || record->isUnmapped()) return;

	int fileIdx = record->getFileIdx();

	//
	// First test whether chr in chrom names match
	//

	bool hasChr = record->hasChrInChromName();
	testType testChrVal = fileHasChrInChromNames(fileIdx);

	if (testChrVal == UNTESTED) {
		_fileHasChrInChromNames[fileIdx] = hasChr ? YES : NO;
	}
	if ((_allFilesHaveChrInChromNames == YES && !hasChr) || (_allFilesHaveChrInChromNames == NO && hasChr)) {
		nameConventionWarning(record, _fileNames[fileIdx], " has inconsistent naming convention for record:\n");
	}

	if (_allFilesHaveChrInChromNames == UNTESTED) {
		_allFilesHaveChrInChromNames = hasChr ? YES : NO;
	}


	//
	// Now test whether leading zero in chrom names match
	//


	bool zeroVal = record->hasLeadingZeroInChromName(hasChr);
	testChrVal = fileHasLeadingZeroInChromNames(fileIdx);
	if (testChrVal == UNTESTED) {
		_fileHasLeadingZeroInChromNames[fileIdx] = zeroVal ? YES : NO;
	}
	if ((_allFileHaveLeadingZeroInChromNames == YES && !zeroVal) || (_allFileHaveLeadingZeroInChromNames == NO && zeroVal)) {
		nameConventionWarning(record, _fileNames[fileIdx], " has a record where naming convention (leading zero) is inconsistent with other files:\n");
	}

	if (_allFileHaveLeadingZeroInChromNames == UNTESTED) {
		_allFileHaveLeadingZeroInChromNames = zeroVal ? YES : NO;
	}
}

ContextBase::testType ContextBase::fileHasChrInChromNames(int fileIdx) {
	conventionType::iterator iter = _fileHasChrInChromNames.find(fileIdx);
	if (iter == _fileHasChrInChromNames.end()) {
		return UNTESTED;
	}
	return iter->second;
}

ContextBase::testType ContextBase::fileHasLeadingZeroInChromNames(int fileIdx) {
	conventionType::iterator iter = _fileHasLeadingZeroInChromNames.find(fileIdx);
	if (iter == _fileHasLeadingZeroInChromNames.end()) {
		return UNTESTED;
	}
	return iter->second;
}



void ContextBase::warn(const Record *record, const string str1, const string str2, const string str3) {
	string msg;
	setErrorMsg(msg, true, record, str1, str2, str3);
	cerr << msg << endl;
}

void ContextBase::die(const Record *record, const string str1, const string str2, const string str3) {
	string msg;
	setErrorMsg(msg, false, record, str1, str2, str3);
	cerr << msg << endl;
	exit(1);
}

void ContextBase::setErrorMsg(string &msg, bool onlyWarn, const Record * record, string str1, const string str2, const string str3) {
	if (onlyWarn) {
		msg = "\n***** WARNING: ";
	} else {
		msg = "\n***** ERROR: ";
	}
	msg.append(str1);
	msg.append(str2);
	msg.append(str3);
	msg.append(" Exiting...\n");
	if (record != NULL) {
		record->print(msg);
	}
}

void ContextBase::nameConventionWarning(const Record *record, const string &filename, const string &message)
{
	 _nameConventionWarningMsg = "***** WARNING: File ";
	 _nameConventionWarningMsg.append(filename);
	 _nameConventionWarningMsg.append(message);
	 record->print(_nameConventionWarningMsg);
	 _nameConventionWarningMsg.append("\n");
	 _nameConventionWarningTripped = true;

	 cerr << _nameConventionWarningMsg << endl;
}

bool ContextBase::strandedToolSupported() {
	//Test that all files have strands. Should be called if any tool
	// invoked with sameStrand / diffStrand option.
	for (int i=0; i < getNumInputFiles(); i++) {

		// make sure file has strand.
		if (!getFile(i)->recordsHaveStrand()) {
			_errorMsg = "\n***** ERROR: stranded ";
			_errorMsg += _origProgramName;
			_errorMsg += " requested, but input file ";
			_errorMsg  += getInputFileName(i);
			_errorMsg  += " does not have strands. *****";
			return false;
		}
		//make sure file is not VCF.
		if (getFile(i)->getFileType() == FileRecordTypeChecker::VCF_FILE_TYPE) {
			_errorMsg = "\n***** ERROR: stranded ";
			_errorMsg += _origProgramName;
			_errorMsg += " not supported for VCF file ";
			_errorMsg += getInputFileName(i);
			_errorMsg += ". *****";
			return false;
		}
	}
	return true;
}

