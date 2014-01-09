/*
 * ContextBase.h
 *
 *  Created on: Feb 11, 2013
 *      Author: nek3d
 */

#ifndef CONTEXTBASE_H_
#define CONTEXTBASE_H_


// The Context class handles the settings for an operation,
// such as merge, intersect, jaccard, etc.
//
// Settings include the input and output parameters, such as input
// files, file types (if explicitly provided), genome files,
// run options, output format, etc.

#include <cstdlib>
#include "version.h"
#include "BedtoolsTypes.h"
#include "FileRecordTypeChecker.h"
#include "NewGenomeFile.h"
#include "api/BamReader.h"
#include "api/BamAux.h"


class ContextBase {
public:
	ContextBase();
	virtual ~ContextBase();

	typedef FileRecordTypeChecker::FILE_TYPE ContextFileType;
	typedef FileRecordTypeChecker::RECORD_TYPE ContextRecordType;

	typedef enum {UNSPECIFIED_PROGRAM, INTERSECT, WINDOW, CLOSEST, COVERAGE, MAP, GENOMECOV, MERGE, CLUSTER,
		COMPLEMENT, SUBTRACT, SLOP, FLANK, SORT, RANDOM, SAMPLE, SHUFFLE, ANNOTATE, MULTIINTER, UNIONBEDG, PAIRTOBED,
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
	//HERE ARE SOME SIMPLER VERSIONS OF THE ABOVE FOR APPS THAT HAVE ONLY ONE INPUT FILE
	const QuickString &getInputFileName() const { return _inputFiles[0]._fileName; }
	ContextFileType getInputFileType() const { return _inputFiles[0]._fileType; }
	void setInputFileType(ContextFileType fileType) { _inputFiles[0]._fileType = fileType; }
	ContextRecordType getInputRecordType() const { return _inputFiles[0]._recordType; }
	void setInputRecordType(ContextRecordType recordType) { _inputFiles[0]._recordType = recordType; }
	int getInputFileIdx() const { return 0; }

	virtual bool determineOutputType();

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

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);

	 //isValidState checks that parameters to context are in an acceptable state.
	// If not, the error msg string will be set with a reason why it failed.
	virtual bool isValidState();
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

    bool getUseBufferedOutput() const { return _useBufferedOutput; }
    void setUseBufferedOutput(bool val) { _useBufferedOutput = val; }

    virtual bool getSortedInput() const {return _sortedInput; }
    virtual void setSortedInput(bool val) { _sortedInput = val; }

    virtual bool getPrintHeader() const {return _printHeader; }
    virtual void setPrintHeader(bool val) { _printHeader = val; }

    virtual bool getPrintable() const { return _printable; }
    virtual void setPrintable(bool val) { _printable = val; }

    virtual bool getUseFullBamTags() const { return _useFullBamTags; }
    virtual void setUseFullBamTags(bool val) { _useFullBamTags = val; }

	//
	// MERGE METHODS
	//
    virtual bool getReportCount() const { return _reportCount; }
    virtual void setReportCount(bool val) { _reportCount = val; }

    virtual int getMaxDistance() const { return _maxDistance; }
    virtual void setMaxDistance(int distance) { _maxDistance = distance; }

    virtual bool getReportNames() const { return _reportNames; }
    virtual void setReportNames(bool val) { _reportNames = val; }

    virtual bool getReportScores() const { return _reportScores; }
    virtual void setReportScores(bool val) { _reportScores = val; }

    virtual const QuickString &getScoreOp() const { return _scoreOp; }
    virtual void setScoreOp(const QuickString &op) { _scoreOp = op; }


	// METHODS FOR PROGRAMS WITH USER_SPECIFIED NUMBER
	// OF OUTPUT RECORDS.
    virtual int getNumOutputRecords() const { return _numOutputRecords; }
    virtual void setNumOutputRecords(int val) { _numOutputRecords = val; }

	//SEED OPS FOR APPS WITH RANDOMNESS
	//When a seed has been specified on the command line, srand
	//will have been already been called during command line parsing.
	//For reference, srand is the C function that initializes the
	//psuedo-random number generator. If a seed has not been
	//specifified, the user may call the
    virtual bool hasConstantSeed() const { return _hasConstantSeed; }
    virtual int getConstantSeed() const { return _seed; }
    virtual int getUnspecifiedSeed();

    //method used by other classes to distinguish ContextBase from
    //ContextIntersect and other derived classes that use intersection
    //methods.
    virtual bool hasIntersectMethods() const { return false; }

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
	int _skipFirstArgs;
	bool _showHelp;
    bool _obeySplits;
    bool _uncompressedBam;
    bool _useBufferedOutput;

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
    string _columnOperation;
    int _column;
	bool _reportCount;
	int _maxDistance;
	bool _reportNames;
	bool _reportScores;
	QuickString _scoreOp;
	set<QuickString> _validScoreOps;

	int _numOutputRecords;

	bool _hasConstantSeed;
	int _seed;
	bool _forwardOnly;
	bool _reverseOnly;

	void markUsed(int i) { _argsProcessed[i] = true; }
	bool isUsed(int i) const { return _argsProcessed[i]; }
	bool cmdArgsValid();

	//set cmd line params and counter, i, as members so code
	//is more readable (as opposed to passing all 3 everywhere).
	int _argc;
	char **_argv;
	int _i;

	virtual bool handle_bed();
	virtual bool handle_fbam();
	virtual bool handle_g();
	virtual bool handle_h();
	virtual bool handle_header();
	virtual bool handle_i();
	virtual bool handle_n();
	virtual bool handle_nobuf();
	virtual bool handle_seed();
	virtual bool handle_split();
	virtual bool handle_sorted();
	virtual bool handle_ubam();
};

#endif /* CONTEXTBASE_H_ */
