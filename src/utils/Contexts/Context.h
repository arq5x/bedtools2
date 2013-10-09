/*
 * ContextBase.h
 *
 *  Created on: Feb 11, 2013
 *      Author: nek3d
 */

#ifndef CONTEXTBASE_H_
#define CONTEXTBASE_H_


// This is an abstract base class for all context objects.
//
// Context classes handle the settings for an operation,
// such as merge, intersect, jaccard, etc.
//
// Settings include the input and output parameters, such as input
// files, file types (if explicitly provided), genome files,
// run options, output format, etc.

#include "BedtoolsTypes.h"
#include "FileRecordTypeChecker.h"
#include "NewGenomeFile.h"
#include "api/BamReader.h"
#include "api/BamAux.h"


class Context {
public:
	Context();
	~Context();

	typedef FileRecordTypeChecker::FILE_TYPE ContextFileType;
	typedef FileRecordTypeChecker::RECORD_TYPE ContextRecordType;

	typedef enum {UNSPECIFIED_PROGRAM, INTERSECT, WINDOW, CLOSEST, COVERAGE, MAP, GENOMECOV, MERGE, CLUSTER,
		COMPLEMENT, SUBTRACT, SLOP, FLANK, SORT, RANDOM, SHUFFLE, ANNOTATE, MULTIINTER, UNIONBEDG, PAIRTOBED,
		PAIRTOPAIR,BAMTOBED, BEDTOBAM, BEDTOFASTQ, BEDPETOBAM, BED12TOBED6, GETFASTA, MASKFASTA, NUC,
		MULTICOV, TAG, JACCARD, OVERLAP, IGV, LINKS,MAKEWINDOWS, GROUPBY, EXPAND } PROGRAM_TYPE;

	PROGRAM_TYPE getProgram() const { return _program; }
	void setProgram(PROGRAM_TYPE program) { _program = program; }

	void addInputFile(const QuickString &inputFile,
			ContextFileType explicitFileType = FileRecordTypeChecker::UNKNOWN_FILE_TYPE,
			ContextRecordType explicitRecordType = FileRecordTypeChecker::UNKNOWN_RECORD_TYPE) {
		_inputFiles.push_back(FileEntryType(inputFile, explicitFileType, explicitRecordType));
	}

	int getNumInputFiles() const { return _inputFiles.size(); }
	const QuickString &getInputFileName(int fileNum) const { return _inputFiles[fileNum]._fileName; }
	ContextFileType getInputFileType(int fileNum) const { return _inputFiles[fileNum]._fileType; }
	void setInputFileType(int fileNum, ContextFileType fileType) { _inputFiles[fileNum]._fileType = fileType; }
	ContextRecordType getInputRecordType(int fileNum) const { return _inputFiles[fileNum]._recordType; }
	void setInputRecordType(int fileNum, ContextRecordType recordType) { _inputFiles[fileNum]._recordType = recordType; }


	bool determineOutputType();

	const QuickString &getHeader(int fileIdx) { return _headers[fileIdx]; }
	void setHeader(int fileIdx, const QuickString &header) { _headers[fileIdx] = header; }
	const BamTools::RefVector &getReferences(int fileIdx)  { return _references[fileIdx]; }
	void setReferences(int fileIdx, const BamTools::RefVector &refs) { _references[fileIdx] = refs; }
	int getBamHeaderAndRefIdx(); //return idx of 1st query that is BAM. If none, first DB that is BAM.

	bool getUseMergedIntervals() const { return _useMergedIntervals; }
	void setUseMergedIntervals(bool val) { _useMergedIntervals = val; }

	void openGenomeFile(const QuickString &genomeFilename);
	void openGenomeFile(const BamTools::RefVector &refVector);
	bool hasGenomeFile() const { return _genomeFile != NULL; }
	NewGenomeFile *getGenomeFile() const { return _genomeFile; }

	void setOutputFileType(ContextFileType fileType) { _outputFileType = fileType; }
	ContextFileType getOutputFileType() const { return _outputFileType; }

	bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);

	 //isValidState checks that parameters to context are in an acceptable state.
	// If not, the error msg string will be set with a reason why it failed.
	bool isValidState();
	bool getShowHelp() const { return _showHelp; }
	void setShowHelp(bool val) { _showHelp = val; }

	const QuickString &getErrorMsg() const { return _errorMsg; }
	void setErrorMessage(const QuickString &errorMsg) { _errorMsg = errorMsg; }

	//split handling.
	bool getObeySplits() const {return _obeySplits; }
	void setObeySplits(bool val) { _obeySplits = val; }

	//Decide whether output format is BED or BAM.
	//Default is BAM if any input files are BAM.
    bool getExplicitBedOutput() const { return _explicitBedOutput; }
    void setExplicitBedOutput(bool val) { _explicitBedOutput = val; }

    bool getUncompressedBam() const { return _uncompressedBam; }
    void setUncompressedBam(bool val) { _uncompressedBam = val; }
	//
	// INTERSECT METOHDS
	//
	//NOTE: Query and database files will only be marked as such by either the
	//parseCmdArgs method, or by explicitly setting them.
    int getQueryFileIdx() const { return _queryFileIdx; }
	void setQueryFileIdx(int idx) { _queryFileIdx = idx; }
	int getDatabaseFileIdx() const { return _databaseFileIdx; }
	void setDatabaseFileIdx(int idx) { _databaseFileIdx = idx; }
	const QuickString &getQueryFileName() const { return _inputFiles[_queryFileIdx]._fileName; }
	const QuickString &getDatabaseFileName() const { return _inputFiles[_databaseFileIdx]._fileName; }
	ContextFileType getQueryFileType() const { return _inputFiles[_queryFileIdx]._fileType; }
	ContextFileType getDatabaseFileType() const { return _inputFiles[_databaseFileIdx]._fileType; }
	ContextRecordType getQueryRecordType() const { return _inputFiles[_queryFileIdx]._recordType; }
	ContextRecordType getDatabaseRecordType() const { return _inputFiles[_databaseFileIdx]._recordType; }
	int getMaxNumDatabaseFields() const { return _maxNumDatabaseFields; }
	void setMaxNumDatabaseFields(int val) { _maxNumDatabaseFields = val; }

	bool getAnyHit() const {return _anyHit; }
	void setAnyHit(bool val) { _anyHit = val; }

	bool getNoHit() const {return _noHit; }
	void setNoHit(bool val) { _noHit = val; }

	bool getLeftJoin() const {return _leftJoin; }
	void setLeftJoin(bool val) { _leftJoin = val; }

	bool getWriteA() const {return _writeA; }
	void setWriteA(bool val) { _writeA = val; }

	bool getWriteB() const {return _writeB; }
	void setWriteB(bool val) { _writeB = val; }

	bool getWriteCount() const {return _writeCount; }
	void setWriteCount(bool val) { _writeCount = val; }

	bool getWriteOverlap() const {return _writeOverlap; }
	void setWriteOverlap(bool val) { _writeOverlap = val; }

	bool getWriteAllOverlap() const {return _writeAllOverlap; }
	void setWriteAllOverlap(bool val) { _writeAllOverlap = val; }

	bool getHaveFraction() const {return _haveFraction; }
	void setHaveFraction(bool val) { _haveFraction = val; }

	float getOverlapFraction() const { return _overlapFraction; }
	void setOverlapFraction(float fraction) { _overlapFraction = fraction; }

	bool getReciprocal() const {return _reciprocal; }
	void setReciprocal(bool val) { _reciprocal = val; }

	bool getSameStrand() const {return _sameStrand; }
	void setSameStrand(bool val) { _sameStrand = val; }

	bool getDiffStrand() const {return _diffStrand; }
	void setDiffStrand(bool val) { _diffStrand = val; }

	bool getSortedInput() const {return _sortedInput; }
	void setSortedInput(bool val) { _sortedInput = val; }

	bool getPrintHeader() const {return _printHeader; }
	void setPrintHeader(bool val) { _printHeader = val; }

	bool getPrintable() const { return _printable; }
	void setPrintable(bool val) { _printable = val; }

	bool getUseFullBamTags() const { return _useFullBamTags; }
	void setUseFullBamTags(bool val) { _useFullBamTags = val; }

	//
	// MERGE METHODS
	//
	bool getReportCount() const { return _reportCount; }
	void setReportCount(bool val) { _reportCount = val; }

	int getMaxDistance() const { return _maxDistance; }
	void setMaxDistance(int distance) { _maxDistance = distance; }

	bool getReportNames() const { return _reportNames; }
	void setReportNames(bool val) { _reportNames = val; }

	bool getReportScores() const { return _reportScores; }
	void setReportScores(bool val) { _reportScores = val; }

	const QuickString &getScoreOp() const { return _scoreOp; }
	void setScoreOp(const QuickString &op) { _scoreOp = op; }


protected:
	PROGRAM_TYPE _program;

	class FileEntryType {
	public:
		FileEntryType(const QuickString &filename, ContextFileType fileType, ContextRecordType recordType)
		: _fileName(filename), _fileType(fileType), _recordType(recordType) {}
		QuickString _fileName;
		ContextFileType _fileType;
		ContextRecordType _recordType;
	};
	vector<FileEntryType> _inputFiles;
	map<QuickString, PROGRAM_TYPE> _programNames;

	bool _useMergedIntervals;
	NewGenomeFile *_genomeFile;

	ContextFileType _outputFileType;
	bool _outputTypeDetermined;

	map<int, QuickString> _headers;
	map<int, BamTools::RefVector> _references;

	QuickString _errorMsg;
	vector<bool> _argsProcessed; //used for processing cmdLine args.
	int _argc;
	char **_argv;
	int _skipFirstArgs;
	bool _showHelp;
    bool _obeySplits;
    bool _uncompressedBam;

	bool _anyHit;
    bool _noHit;
    bool _writeA;
    bool _writeB;
    bool _leftJoin;
    bool _writeCount;
    bool _writeOverlap;
    bool _writeAllOverlap;
    bool _haveFraction;
    float _overlapFraction;
    bool _reciprocal;
    bool _sameStrand;
    bool _diffStrand;
    bool _sortedInput;
    bool _printHeader;
    bool _printable;
    bool _explicitBedOutput;
    int _queryFileIdx;
    int _databaseFileIdx;
    int _bamHeaderAndRefIdx;
    int _maxNumDatabaseFields;
    bool _useFullBamTags;

	bool _reportCount;
	int _maxDistance;
	bool _reportNames;
	bool _reportScores;
	QuickString _scoreOp;
	set<QuickString> _validScoreOps;

	void markUsed(int i) { _argsProcessed[i] = true; }
	bool isUsed(int i) const { return _argsProcessed[i]; }
	bool cmdArgsValid();

};

#endif /* CONTEXT_H_ */
