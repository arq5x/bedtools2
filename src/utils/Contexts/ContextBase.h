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
#include "FileRecordMergeMgr.h"
#include "NewGenomeFile.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "KeyListOps.h"


class ContextBase {
public:
	ContextBase();
	virtual ~ContextBase();

	typedef FileRecordTypeChecker::FILE_TYPE ContextFileType;
	typedef FileRecordTypeChecker::RECORD_TYPE ContextRecordType;

	typedef enum {UNSPECIFIED_PROGRAM, INTERSECT, WINDOW, CLOSEST, COVERAGE, MAP, GENOMECOV, MERGE, CLUSTER,
		COMPLEMENT, SUBTRACT, SLOP, FLANK, SORT, RANDOM, SAMPLE, SHUFFLE, ANNOTATE, MULTIINTER, UNIONBEDG, PAIRTOBED,
		PAIRTOPAIR,BAMTOBED, BEDTOBAM, BEDTOFASTQ, BEDPETOBAM, BED12TOBED6, GETFASTA, MASKFASTA, NUC,
		MULTICOV, TAG, JACCARD, OVERLAP, IGV, LINKS,MAKEWINDOWS, GROUPBY, EXPAND, SPACING, FISHER, GROUP_BY} PROGRAM_TYPE;

	PROGRAM_TYPE getProgram() const { return _program; }
	FileRecordMgr *getFile(int fileIdx) { return _files[fileIdx]; }
	void setProgram(PROGRAM_TYPE program) { _program = program; }

	void addInputFile(const string &inputFile) { _fileNames.push_back(inputFile); }

	int getNumInputFiles() const { return _fileNames.size(); }
	const string &getInputFileName(int fileNum) const { return _fileNames[fileNum]; }
	ContextFileType getInputFileType(int fileNum) const { return _files[fileNum]->getFileType(); }
	ContextRecordType getInputRecordType(int fileNum) const { return _files[fileNum]->getRecordType(); }

	virtual bool determineOutputType();

	const string &getHeader(int fileIdx) { return _files[fileIdx]->getHeader(); }
	const BamTools::RefVector &getBamReferences(int fileIdx)  { return _files[fileIdx]->getBamReferences(); }
	refs_t* getCramRefs(int fileIdx) { return _files[fileIdx]->getCramRefs(); }
	int getBamHeaderAndRefIdx(); //return idx of 1st query that is BAM. If none, first DB that is BAM.

	bool getUseMergedIntervals() const { return _useMergedIntervals; }
	void setUseMergedIntervals(bool val) { _useMergedIntervals = val; }
	FileRecordMergeMgr::WANTED_STRAND_TYPE getDesiredStrand() const { return _desiredStrand; }

	void openGenomeFile(const string &genomeFilename);
	void openGenomeFile(const BamTools::RefVector &refVector);
	bool hasGenomeFile() const { return _genomeFile != NULL; }
	NewGenomeFile *getGenomeFile() const { return _genomeFile; }

	void setOutputFileType(ContextFileType fileType) { _outputFileType = fileType; }
	ContextFileType getOutputFileType() const { return _outputFileType; }

	virtual bool testCmdArgs(int argc, char **argv);
	virtual bool errorEncountered();

	 //isValidState checks that parameters to context are in an acceptable state.
	// If not, the error msg string will be set with a reason why it failed.
	virtual bool isValidState();
	bool getShowHelp() const { return _showHelp; }
	void setShowHelp(bool val) { _showHelp = val; }

	const string &getErrorMsg() const { return _errorMsg; }
	void setErrorMessage(const string &errorMsg) { _errorMsg = errorMsg; }

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

    virtual bool getSortOutput() const {return _sortOutput; }
    virtual void setSortOutput(bool val) { _sortOutput = val; }

    virtual bool getUseDBnameTags() const { return _reportDBnameTags; }
    virtual void setUseDBnameTags(bool val) { _reportDBnameTags = val; }

    virtual bool getUseDBfileNames() const { return _reportDBfileNames; }
    virtual void setUseDBfileNames(bool val) { _reportDBfileNames = val; }

    virtual bool getPrintHeader() const {return _printHeader; }
    virtual void setPrintHeader(bool val) { _printHeader = val; }

    virtual bool getPrintable() const { return _printable; }
    virtual void setPrintable(bool val) { _printable = val; }

    virtual bool getUseFullBamTags() const { return _useFullBamTags; }
    virtual void setUseFullBamTags(bool val) { _useFullBamTags = val; }

    virtual bool getNameCheckDisabled() const { return _nameCheckDisabled; }
    virtual void setNameCheckDisabled(bool val) { _nameCheckDisabled = val; }

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

    // determine whether column operations like those used in map
    // are available.
    void setColumnOpsMethods(bool val);
    virtual bool hasColumnOpsMethods() const { return _hasColumnOpsMethods; }
    const string &getColumnOpsVal(RecordKeyVector &keyList) const;
    //methods applicable only to column operations.
    int getReportPrecision() const { return _reportPrecision; }


    void testNameConventions(const Record *);
    BlockMgr *getSplitBlockInfo() { return _splitBlockInfo; }

    string getErrorMessages() const { return _errorMsg; }

	bool isCram() const { return _isCram; }

protected:
	PROGRAM_TYPE _program;

	vector<string> _fileNames;
	vector<FileRecordMgr *> _files;
	bool _allFilesOpened;
	map<string, PROGRAM_TYPE> _programNames;
	string _origProgramName;

	NewGenomeFile *_genomeFile;

	ContextFileType _outputFileType;
	bool _outputTypeDetermined;

	map<int, string> _headers;
	map<int, BamTools::RefVector> _references;

	string _errorMsg;
	vector<bool> _argsProcessed; //used for processing cmdLine args.
	int _skipFirstArgs;
	bool _showHelp;
    bool _obeySplits;
    bool _uncompressedBam;
    bool _useBufferedOutput;
    size_t _ioBufSize;

	bool _anyHit;
    bool _noHit;
    bool _writeA;
    bool _writeB;
    bool _leftJoin;
    bool _writeCount;
    bool _writeCountsPerDatabase;
    bool _writeOverlap;
    bool _writeAllOverlap;
    bool _haveFractionA;
    bool _haveFractionB;
    float _overlapFractionA;
    float _overlapFractionB;
    bool _reciprocalFraction;
    bool _eitherFraction;
    bool _sameStrand;
    bool _diffStrand;
    bool _sortedInput;
    bool _sortOutput;
    bool _reportDBnameTags;
    bool _reportDBfileNames;
    bool _printHeader;
    bool _printable;
    bool _explicitBedOutput;
    bool _runToQueryEnd;
    int _queryFileIdx;
    vector<int> _dbFileIdxs;
    vector<string> _dbNameTags;
    map<int, int> _fileIdsToDbIdxs;
    int _bamHeaderAndRefIdx;
    int _maxNumDatabaseFields;
    bool _useFullBamTags;
	bool _isCram;   // Used when a "BAM" type which is actually a CRAM

	int _numOutputRecords;

	bool _hasConstantSeed;
	int _seed;
	bool _forwardOnly;
	bool _reverseOnly;
	bool _nameCheckDisabled;

	//Members for column operations
	bool _hasColumnOpsMethods;
	KeyListOps *_keyListOps;
	string _nullStr; //placeholder return value when col ops aren't valid.

	//Members for merged records
	FileRecordMergeMgr::WANTED_STRAND_TYPE _desiredStrand;
	int _maxDistance;
	bool _useMergedIntervals;

	int _reportPrecision; //used in fields reported from numeric ops from map and merge.


	void markUsed(int i) { _argsProcessed[i] = true; }
	bool isUsed(int i) const { return _argsProcessed[i]; }
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	bool cmdArgsValid();
	virtual bool openFiles();
	virtual FileRecordMgr *getNewFRM(const string &filename, int fileIdx);

	//set cmd line params and counter, i, as members so code
	//is more readable (as opposed to passing all 3 everywhere).
	int _argc;
	char **_argv;
	int _i;

	//track whether each file has the letters chr in their chrom names.
	//this is needed for enforcing consistent naming conventions across
	//multiple files, and auto-detecting errors when they're not followed.
	typedef enum { YES, NO, UNTESTED } testType;
	typedef map<int, testType> conventionType;
	conventionType _fileHasChrInChromNames;
	// as above, but check whether first digit to appear in
	// a chrom name is a zero.
	conventionType _fileHasLeadingZeroInChromNames;

	static const int MIN_ALLOWED_BUF_SIZE = 8;
	BlockMgr *_splitBlockInfo;

    testType _allFilesHaveChrInChromNames;
    testType _allFileHaveLeadingZeroInChromNames;
    bool _noEnforceCoordSort;
	bool _inheader;


    virtual bool handle_bed();
	virtual bool handle_fbam();
	virtual bool handle_g();
	virtual bool handle_h();
	virtual bool handle_header();
	virtual bool handle_i();
	virtual bool handle_n();
	virtual bool handle_nobuf();
	virtual bool handle_iobuf();

	virtual bool handle_seed();
	virtual bool handle_split();
	virtual bool handle_sorted();
	virtual bool handle_ubam();

	virtual bool handle_c();
	virtual bool handle_o();
	virtual bool handle_null();
	virtual bool handle_delim();
	virtual bool handle_sortout();
	virtual bool handle_nonamecheck();
	bool handle_prec();
	bool parseIoBufSize(string bufStr);

    testType fileHasChrInChromNames(int fileIdx);
    testType fileHasLeadingZeroInChromNames(int fileIdx);

    void setNoEnforceCoordSort(bool val) { _noEnforceCoordSort = val; }
    //Warning messages.
   bool _nameConventionWarningTripped;
   string _nameConventionWarningMsg;
   void nameConventionWarning(const Record *record, const string &filename, const string &message);

    //give warning but continue.
    void warn(const Record *, const string str1, const string str2 = "", const string str3 = "");
    // Give error and exit.
    void die(const Record *, const string str1, const string str2 = "", const string str3 = "");

    //private error handler
    void setErrorMsg(string &msg, bool onlyWarn, const Record * record, string str1, const string str2, const string str3);

    bool strandedToolSupported();

};

#endif /* CONTEXTBASE_H_ */
