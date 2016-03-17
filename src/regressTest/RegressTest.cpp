#include "RegressTest.h"
#include <cstdlib>
#include "SingleLineDelimTextFileReader.h"
#include "BufferedStreamMgr.h"
#include "lineFileUtilities.h"
#include "ParseTools.h"
#include <sys/stat.h>
#include <ctime>

const string RegressTest::_hardOptsCmd = "HARD_OPTIONS";
const string RegressTest::_correctCmd = "CORRECT";
const string RegressTest::_performCmd = "PERFORM";
const string RegressTest::_randomCmd = "RANDOM";

const string RegressTest::_space = " ";
const string RegressTest::_redirect = " > ";
const string RegressTest::_devNull = " /dev/null ";
//const string RegressTest::_bedOpsCmd = "/home/nek3d/testWorkspace/bedops/bin/bedmap --echo  --echo-map  --bp-ovr 1 ";
const string RegressTest::_bedOpsCmd = "/home/nek3d/testWorkspace/bedops/bin/bedops --intersect ";


RegressTest::RegressTest()
:   _filesPerRun(0),
    _generatedFileNumber(1),
    _isOldProgBedops(false)
{
	_correctFiles = new fileListType();
	_performFiles = new fileListType();
}

RegressTest::~RegressTest()
{

	delete _correctFiles;
	_correctFiles = NULL;

	delete _performFiles;
	_performFiles = NULL;
}

bool RegressTest::init(int argc, char **argv)
{
	//TBD: Perhaps allow option for specifying report file.
	_fpReportFile = stdout;

	if (!parseParams(argc, argv)) {
		return false;
	}

	//make banner title for report
	fprintf(_fpReportFile, "\n\n***********************************************************\n\n");
	fprintf(_fpReportFile, "REGRESSION TEST FOR BEDTOOLS SUB-PROGRAM: %s\n", _subProgram.c_str());
	fprintf(_fpReportFile, "\n***********************************************************\n\n");


	//setup output directory for testing's temporary files and generated data files.
	_tmpDirname = "tempTesting_";
	time_t currTime = time(NULL);
	string timeStr = ctime(&currTime); //timeStr will equal Www Mmm dd hh:mm:ss yyyy followed by newline and null char.
	//chomp the newLine
	timeStr.erase(timeStr.size()-1);
	//adjust timeStr to change spaces to underscores.
	for (int i=0; i < (int)timeStr.size(); i++) {
		if (isspace(timeStr[i])) {
			timeStr[i] = '_';
		} else if (timeStr[i] == ':') {
			timeStr[i] = '-';
		}
	}
	_tmpDirname += timeStr;
	int mkdirRetval = mkdir(_tmpDirname.c_str(), S_IRWXU | S_IRWXG | S_IRWXO ); //mkdir directory with all permissions allowed.
	if (mkdirRetval != 0) {
		fprintf(stderr, "Error: Unable to create temporary output directory %s.\n", _tmpDirname.c_str());
		return false;
	}
	_tmpDirname += "/";

	_memoryLogfilename = _tmpDirname + "memoryLog.txt";

	_userName = getenv("USER");

	if ((int)_filePrecessorOptions.size() != _filesPerRun) {
		fprintf(stderr, "Error: Sub-program specific filesPerRun number must match number of precessor options.\n");
		return false;
	}

	if (!config()) {
		return false;
	}

	if ((int)_correctFiles->size() % _filesPerRun != 0) {
		fprintf(stderr, "Error: number of correctness files is not evenly divisible by number of files per run.\n");
		return false;
	}

	if ((int)_performFiles->size() % _filesPerRun != 0) {
		fprintf(stderr, "Error: number of performance files is not evenly divisible by number of files per run.\n");
		return false;
	}

	return true;
}

bool RegressTest::runTests() {
	echoOptions();

	//All set with set-up. Let's work some mojo.
	if (!performTests(true)) {
		fprintf(stderr, "Error: Failure in correctness tests.\n");
		return false;
	}

	if (!performTests(false)) {
		fprintf(stderr, "Error: Failure in performance tests.\n");
		return false;
	}

	return true;

}

bool RegressTest::parseParams(int argc, char **argv)
{
	//usage: bedtools devtest sub-prog targetVersion configFile [optionsToTest]
	if (argc < 5) {
		usage();
		return false;
	}
	_newVersion = argv[0];
	_subProgram = argv[2];

	_oldVersion = argv[3];
	if (_oldVersion.find("bedops") != string::npos) {
		_isOldProgBedops = true;
	}


	_configFilename = argv[4];

	//special: add a blank option to the softOptions, so that performTests will compare
	//runs with each soft opt to a run with no soft ops.
	_softOptions.push_back("");
	for (int i=5; i < argc; i++ ) {
		_softOptions.push_back(argv[i]);
	}
	return true;
}

void RegressTest::usage() const {
	fprintf(stderr, "Usage: bedtools regresstest sub-prog targetVersion configFile [optionsToTest]\n");
}

void RegressTest::echoOptions() const
{
	//show all command line and config file options and files
	fprintf(_fpReportFile, "\nCONFIGURATION AND OPTIONS ARE:\n\n");
	fprintf(_fpReportFile, "New Version: %s\n", _newVersion.c_str());
	fprintf(_fpReportFile, "Old Version: %s\n\n", _oldVersion.c_str());
	fprintf(_fpReportFile, "hardOptions: %s\n\n", _hardOptions.c_str());

	for (int i=1; i < (int)_softOptions.size(); i++) {
		fprintf(_fpReportFile, "SoftOption %d: %s\n", i, _softOptions[i].c_str());
	}

	fprintf(_fpReportFile, "\nFiles per run: %d\n", _filesPerRun);
	for (int i=0; i < (int)_filePrecessorOptions.size(); i++) {
		fprintf(_fpReportFile, "Precessor option %d: %s\n", i+1, _filePrecessorOptions[i].c_str());
	}

	for (int i=0; i < 2; i++) {
		fileListType *fileList = i == 0 ?  _correctFiles : _performFiles;
		fprintf(_fpReportFile, "\n\n%s TEST FILES:\n", i ==0 ? "CORRECTNESS" : "PERFORMANCE");
		for (fileListType::const_iterator iter = fileList->begin(); iter != fileList->end(); iter++) {

			const string &filename = iter->first;
			const string &desc = iter->second;
			fprintf(_fpReportFile, "\nFilename:  %s\n", filename.c_str());
			fprintf(_fpReportFile, "Description: %s\n", desc.c_str());
		}
	}
}

bool RegressTest::config() {
	//Set hard options, populate correctness and performance file vectors by reading a config file

	_configFile.open(_configFilename.c_str());
	if (!_configFile.good()) {
		cerr << "Error: unable to open config file " << _configFilename << endl;
	}
	vector<string> fields;
	bool parseStatus = false;
	int numLinesToRead = 0;
	int linesRead = 0;
	string sLine;
	while (!_configFile.eof()) {
		sLine.clear();
		getline(_configFile, sLine);
		fields.clear();

		Tokenize(sLine, fields);
		if (fields.size() != 2) {
			continue;
		}

		const string &field1 = fields[0];
		const string &field2 = fields[1];

		if (field1 == _hardOptsCmd) {
			_hardOptions = field2.c_str();
			continue;
		} else if (field1 == _correctCmd || field1 == _performCmd) {
			numLinesToRead = atoi(field2.c_str());
			if (field1 == _correctCmd) {
				parseStatus = parseConfigLines(numLinesToRead, true);
			} else {
				parseStatus = parseConfigLines(numLinesToRead, false);
			}
			if (!parseStatus) {
				fprintf(stderr, "Error: failed to read and parse requested %d lines for %s\n",
						numLinesToRead, field1 == _correctCmd ? "Correctness tests" : "Performance tests");
				_configFile.close();
				return false;
			} else {
				linesRead += numLinesToRead;
			}
		} else {
			fprintf(stderr, "Error: Malformed config file %s\n.\tCheck that num files specified matches num provided.\n", _configFilename.c_str());
			_configFile.close();
			return false;
		}
	}
	if (linesRead == 0) {
		fprintf(stderr, "Error: No file lines read in config file %s\n", _configFilename.c_str());
		_configFile.close();
		return false;
	}

	_configFile.close();
	return true;
}



bool RegressTest::parseConfigLines(int numLinesToRead, bool correctnessFiles)
{
	string description;


	fileListType *fileList = correctnessFiles ? _correctFiles : _performFiles;
	string sLine;
	vector<string> fields;
	for (int i=0; i < numLinesToRead; i++) {

		sLine.clear();
		getline(_configFile, sLine);
		fields.clear();

		Tokenize(sLine, fields);
		if (fields.size() != 2) {
			continue;
		}

		const string &field1 = fields[0];
		const string &field2 = fields[1];

		if (field1 == _randomCmd) {
			string genFilename;
			if (!generateRandomFile(field2.c_str(), genFilename)) {
				fprintf(stderr, "Error: could not generate random file with args: %s\n", field2.c_str());
				return false;
			}
			description = _randomCmd + _space + field2.c_str();
			fileList->push_back(make_pair(genFilename, description));
		} else {
			fileList->push_back(make_pair(field1.c_str(), field2.c_str()));
		}
	}
	return true;
}

bool RegressTest::generateRandomFile(const string & randomArgs, string &filename)
{
	if (_generatedFileNumber == 1) { //first call, print banner to report
		fprintf(_fpReportFile, "\nGENERATING TEST DATA:\n\n");
	}


	string genFilename = _tmpDirname + "generatedFile_";
	string strNum;
	int2str(_generatedFileNumber, strNum);
	_generatedFileNumber++;

	genFilename += strNum;
	genFilename += ".bed";

	//quick hack since v2.18 random currently doesn't support the -sort option in the random command
	string genFileCmd = "~/testWorkspace/pfm3release/bin/bedtools"; //_newVersion;
	genFileCmd += " random ";
	string sortCmd = ""; //" | sort -k1,1 -nk2,2 "; Trouble getting pipes to work. Use -sorted command from now.
	genFileCmd += randomArgs + sortCmd + _redirect + genFilename;

	fprintf(_fpReportFile, "Creating data file %s with random args %s...\n", genFilename.c_str(), randomArgs.c_str());
	if (system(genFileCmd.c_str()) != 0) { //expects successful calls to randomBed to return zero.
		fprintf(_fpReportFile, "FAILED.\n");
		filename.clear();
		return false;
	}
	filename = genFilename;
	fprintf(_fpReportFile, "Done.\n");
	return true;
}

bool RegressTest::performTests(bool isCorrectnessTest)
{
	bool retval = true;
	fileListType *fileList = isCorrectnessTest ? _correctFiles : _performFiles;
	if (fileList->empty()) {
		return true;
	}
	string inFilesCmd;
	string baseCmd;
	string allOptsCmd;
	string finalCmd;
	string oldVerionFinalCmd;
	string newVersionFinalCmd;

	int outFileCounter = 1;
	string counterStr;
	string outputBase = "tmp";
	string outputSuffix = ".txt";
	string outputFile = "/dev/null"; // /dev/null will be used for performance, but over-written for correctness


	string bedOpsFileCmd; //special handling for bedops

	string testType = isCorrectnessTest ? "CORRECTNESS" : "PERFORMANCE";
	int testCounter = 0;
	fprintf(_fpReportFile, "\n\n\n***********************************************************\n\n");
	fprintf(_fpReportFile, "TESTS FOR %s\n", testType.c_str());
	fprintf(_fpReportFile, "\n***********************************************************\n");

	vector< fileListType::const_iterator> currFiles;
	for (fileListType::const_iterator fileIter = fileList->begin(); fileIter != fileList->end(); fileIter += _filesPerRun) {
		currFiles.clear();
		string inFilesCmd = _space;
		for (int i=0; i < _filesPerRun; i++) {
			const string &infileName = (fileIter + i)->first;
			inFilesCmd += _filePrecessorOptions[i] + _space + infileName + _space;
			currFiles.push_back(fileIter + i);
		}
		if (_isOldProgBedops) {
			bedOpsFileCmd = fileIter->first + _space + (fileIter + 1)->first + _space;
		}
		for (int i=0; i < (int)_softOptions.size(); i++ ) {
			const string &softOpt = _softOptions[i];

			//stop to report what test we're going to do.
			testCounter++;
			fprintf(_fpReportFile, "\n%s TEST %d\n\n", testType.c_str(), testCounter);
			fprintf(_fpReportFile, "\tInput Files:\n");
			for (int j=0; j < (int)currFiles.size(); j++) {
				fprintf(_fpReportFile, "\tName: %s\n", currFiles[j]->first.c_str());
				fprintf(_fpReportFile, "\tDesc: %s\n\n", currFiles[j]->second.c_str());
			}
			fprintf(_fpReportFile, "\tSoft option in use: %s\n", softOpt.empty() ? "NONE" : softOpt.c_str());
			fprintf(_fpReportFile, "\n\tCommands to run:\n");
			fileListType cmdAndOutput;
			for (int j=0; j < 2; j++) { //loop through the two versions, old and new
				const string &mainProg = j == 0 ? _oldVersion : _newVersion;
				if (j==1 || !_isOldProgBedops) {
					baseCmd = mainProg + _space + _subProgram + _space + inFilesCmd + _hardOptions + _space + softOpt + _redirect;
				} else {
					baseCmd = /*mainProg + _space + */ _bedOpsCmd + bedOpsFileCmd + _redirect;
				}

				//now we just need tmp output for correctness, or clock and memory footprinting for performance
				if (isCorrectnessTest) {
					int2str(outFileCounter, counterStr);
					outFileCounter++;
					outputFile = _tmpDirname + outputBase + counterStr + outputSuffix;
				}

				finalCmd = baseCmd + outputFile;
				fprintf(_fpReportFile, "\t%s\n", finalCmd.c_str());
				cmdAndOutput.push_back(make_pair(finalCmd, outputFile));
			}
			if (isCorrectnessTest) {
				bool compStatus = executeAndCompareCorrectness(cmdAndOutput);
				retval = retval && compStatus;
			} else {
				bool compStatus = executeAndComparePerformance(cmdAndOutput);
				retval = retval && compStatus;
			}
		}
	}

	return retval;

}

bool RegressTest::executeAndCompareCorrectness(const fileListType &fileList) {
	const string &cmd1 =fileList[0].first;
	const string &output1 = fileList[0].second;
	const string &cmd2 = fileList[1].first;
	const string &output2 = fileList[1].second;

//	printf("\nRun correctness test: %s\n%s\n%s\n%s\n", cmd1.c_str(), output1.c_str(), cmd2.c_str(), output2.c_str());
//	return true;

	int ret1 = system(cmd1.c_str());
	if (ret1 !=0) {
		fprintf(stderr, "\nError: received non-zero exit code %d from old version.\n", ret1);
		return false;
	}

	int ret2 = system(cmd2.c_str());

	if (ret2 !=0) {
		fprintf(stderr, "\nError: received non-zero exit code %d from new version.\n", ret2);
		return false;
	}

	//Here: implement a way to test differences in output. In the future, may wish to actually read the records
	//and determine equivalence, i.e. same value for core fields chrName, start, end, name, score, strand.
	//Right now, just use diff command, ensure output of diff is empty.

	string diffFilename = _tmpDirname + "diffOut.txt";
	string diffCmd = "diff ";
	diffCmd += output1 + _space + output2 + _redirect + diffFilename;

	system(diffCmd.c_str());

	//now need to check for empty diffFile.
    struct stat buf ;
    int i;

    i = stat(diffFilename.c_str(), &buf);
    if (i!=0) {
       fprintf(stderr, "Error: can't get status of diff output file %s\n", diffFilename.c_str()) ;
       return false;
    }

    if (buf.st_size > 0) {
    	fprintf(_fpReportFile, "\n\tFAILED. Output files are different.\n");
    	return true;
    }

    fprintf(_fpReportFile, "\n\tPASSED. Output files are identical.\n");
    return true;
}

bool RegressTest::executeAndComparePerformance(const fileListType &fileList) {
	const string &cmd1 =fileList[0].first;
//	const string &output1 = fileList[0].second;
	const string &cmd2 = fileList[1].first;
//	const string &output2 = fileList[1].second;

//	printf("\nRun performance test: %s\n%s\n%s\n%s\n", cmd1.c_str(), output1.c_str(), cmd2.c_str(), output2.c_str());
//	return true;

	bool cmd1IsBedops = cmd1.find("bedops") != string::npos;
	bool cmd2IsBedops = cmd2.find("bedops") != string::npos;

	char timeBuf[100]; //be sure to initialize before use.
	//run and time old version
	fprintf(_fpReportFile, "\n\tRunning first command...\n");


	startMemoryProfile(cmd1IsBedops);

	time_t oldStartTime = time(NULL);
	int ret1 = system(cmd1.c_str());
	if (ret1 !=0) {
		fprintf(stderr, "\nError: received non-zero exit code %d from old version.\n", ret1);
		return false;
	}
	time_t oldEndTime = time(NULL);
	time_t oldRunTime = oldEndTime - oldStartTime;
	struct tm *oldRunTimeInfo = gmtime(&oldRunTime);
	memset(timeBuf, 0, 100);
	strftime(timeBuf, 100, "%X", oldRunTimeInfo);
	fprintf(_fpReportFile, "\tDone. Elaspsed time for old version was %s.\n", timeBuf);

	endMemoryProfile();

	calcMemoryStats();



	//run and time new version
	fprintf(_fpReportFile, "\n\tRunning second command...\n");

	startMemoryProfile(cmd2IsBedops);

	time_t newStartTime = time(NULL);
	int ret2 = system(cmd2.c_str());
	if (ret2 !=0) {
		fprintf(stderr, "\nError: received non-zero exit code %d from new version.\n", ret2);
		return false;
	}
	time_t newEndTime = time(NULL);
	time_t newRunTime = newEndTime - newStartTime;
	struct tm *newRunTimeInfo = gmtime(&newRunTime);
	memset(timeBuf, 0, 100);
	strftime(timeBuf, 100, "%X", newRunTimeInfo);
	fprintf(_fpReportFile, "\tDone. Elaspsed time for new version was %s.\n", timeBuf);

	endMemoryProfile();

	calcMemoryStats();

	if (newRunTime < oldRunTime) { //new version is faster
		time_t diffTime = oldRunTime - newRunTime;
		struct tm *diffRunTimeInfo = gmtime(&diffTime);
		memset(timeBuf, 0, 100);
		strftime(timeBuf, 100, "%X", diffRunTimeInfo);

		fprintf(_fpReportFile, "\n\tPASSED. New version is faster by %s", timeBuf);

		if (newRunTime > 0) {
			fprintf(_fpReportFile, ", or %7.2fx faster\n", (float)oldRunTime/(float)newRunTime);
		} else {
			fprintf(_fpReportFile, ", or ...WAY FASTER.\n");
		}

	} else if (oldRunTime < newRunTime) { //old version is faster
			time_t diffTime = newRunTime - oldRunTime;
			struct tm *diffRunTimeInfo = gmtime(&diffTime);

			memset(timeBuf, 0, 100);
			strftime(timeBuf, 100, "%X", diffRunTimeInfo);
			fprintf(_fpReportFile, "\n\tFAILED. old version is faster by %s", timeBuf);

			if (oldRunTime > 0) {
				fprintf(_fpReportFile, ", or %7.2fx faster\n", (float)newRunTime/(float)oldRunTime);
			} else {
				fprintf(_fpReportFile, ", or ...WAY FASTER.\n");
			}

	} else { //run times were same.
		fprintf(_fpReportFile, "\n\tMEH. No difference in speed.\n");
	}

	return true;
}


bool RegressTest::startMemoryProfile(bool isBedops)
{
	//kick off ps as a background process to monitor memory usage
	string pidFilename = _tmpDirname + "pidFile.txt";
	string psCmd;
	if (!isBedops) {
		psCmd = "while [ 1 ] ; do ps au | grep bedtools | grep ";
	} else {
		psCmd = "while [ 1 ] ; do ps au | grep bedops | grep ";
	}
	psCmd += _userName;
	psCmd += " | grep -v grep | grep -v regresstest | grep -v \"sh -c\" >> ";
	psCmd += _memoryLogfilename;
	psCmd += " ; sleep .5 ; done & echo $! > ";
	psCmd += _tmpDirname + "pidFile.txt";

	if (system(psCmd.c_str()) != 0) {
		fprintf(stderr, "Error: unable to launch top as background process for memory profiling.\n");
		return false;
	}

	// open and parse the pidFile to get the process id of the memory profiling background process.
	FILE *fp = fopen(pidFilename.c_str(), "r");
	char sLine[4192];
	memset(sLine, 0, 4192);
	fgets(sLine, 4192, fp);
	fclose(fp);

	_currPidOfMemoryLogging = sLine;
	return true;
}

bool RegressTest::endMemoryProfile()
{

	string killCmd = "kill -9 ";
	killCmd += _currPidOfMemoryLogging;
	if (system(killCmd.c_str()) != 0) {
		fprintf(stderr, "Error: failed to kill process id %s\n", _currPidOfMemoryLogging.c_str());
		return false;
	}

	return true;
}

bool RegressTest::calcMemoryStats()
{
	//read and parse the memory log file, calc basic stats: max, mean, median.
	FILE *fp = fopen(_memoryLogfilename.c_str(), "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: Unable to read memory profile log %s.\n", _memoryLogfilename.c_str());
		return false;
	}

	char sLine[4192];
	int numValidLines=0;
	char field1[2048];
	char field2[2048];
	char field3[2048];
	char field4[2048];
	char field5[2048];

	int totalMem= 0;
	int currMem = 0;
	int maxMem = 0;
	while(!feof(fp)) {
		memset(sLine, 0, 4192);
		memset(field1, 0, 2048);
		memset(field2, 0, 2048);
		memset(field3, 0, 2048);
		memset(field4, 0, 2048);
		memset(field5, 0, 2048);

		currMem = 0;
		fgets(sLine, 4192, fp);
		int len = strlen(sLine);
		bool isValidLine = false;
		for (int i=0; i < len; i++) {
			if (!isspace(sLine[i])) {
				isValidLine = true;
				break;
			}
		}
		if (isValidLine) {
			numValidLines++;
			sscanf(sLine, "%s %s %s %s %s", field1, field2, field3, field4, field5);

			//now field 5 has the number of kilobytes the process used at that moment.
			currMem = atoi(field5);
			if (currMem == 0) {
				//atoi failed. something wrong with input
				fprintf(stderr, "Error: bad field where memory usage expected: %s.\n", field5);
				fclose(fp);
				return false;
			}
			maxMem = max(currMem, maxMem);
			totalMem += currMem;
		}
	}

	int avgMem = numValidLines > 0 ? totalMem / numValidLines : -1;

	fclose(fp);
	if (avgMem > 0) {
		fprintf(_fpReportFile, "\tMemory uage: max = %dkb\tavg = %dkb\n", maxMem, avgMem);
	} else {
		fprintf(_fpReportFile, "\tMemory usage: TOO FAST TO PROFILE!\n");
	}

	//cleanup: erase memoryLogfile, as we wish to append to a blank file when next used.
	if (remove(_memoryLogfilename.c_str()) != 0) {
		fprintf(stderr, "Error: couldn't delete old memory logfile %s.\n", _memoryLogfilename.c_str());
		return false;
	}

	return true;
}
