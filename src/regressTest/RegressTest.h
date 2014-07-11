/*
 * RegressTest.h
 *
 *  Created on: Dec 21, 2012
 *      Author: nek3d
 */

#ifndef REGRESSTEST_H_
#define REGRESSTEST_H_

#include <string>
#include <cstdio>
#include <vector>
#include <fstream>

using namespace std;

class SingleLineDelimTextFileReader;
class BufferedStreamMgr;

class RegressTest {
public:
	RegressTest();
	~RegressTest();
	bool init(int argc, char **argv);

	void setFilesPerRun(int numFiles) { _filesPerRun = numFiles; }
	void addFilePrecessorOption(const char *option) { _filePrecessorOptions.push_back(option); }
	bool runTests();



protected:
	string _configFilename;
	string _newVersion;
	string _subProgram;
	string _oldVersion;

	string _hardOptions;
	vector<string> _softOptions;

	int _filesPerRun;
	vector<string> _filePrecessorOptions;

	FILE *_fpReportFile;
	ifstream _configFile;

	typedef vector<pair<string, string> > fileListType;
	fileListType *_correctFiles; //list of files for correctness tests. First in pair is file name, second is description.
	fileListType *_performFiles; //list of files for performance tests. First in pair is file name, second is description.

	//config file key words
	static const string _hardOptsCmd;
	static const string _correctCmd;
	static const string _performCmd;
	static const string _randomCmd;

	//useful strings to have for building command strings
	static const string _space;
	static const string _redirect;
	static const string _devNull;
	static const string _bedOpsCmd;

	int _generatedFileNumber; // a "tag" to give generatedFiles.

	string _tmpDirname;
	string _memoryLogfilename;
	string _currPidOfMemoryLogging;
	string _userName;

	//Special: Since RegressTest was originally designed only for testing bedtools against prior versions
	// of itself, some hacking is needed to make it play with bedops.
	bool _isOldProgBedops;


	bool parseParams(int argc, char **argv);
	void usage() const;
	void echoOptions() const;
	bool config();
	bool parseConfigLines(int numLinesToRead, bool correctnessFiles);
	bool generateRandomFile(const string & randomArgs, string &filename);

	bool performTests(bool correctness); //pass true for correctness, false for performance.
	bool executeAndCompareCorrectness(const fileListType &cmdAndOutput);
	bool executeAndComparePerformance(const fileListType &cmdAndOutput);

	bool startMemoryProfile(bool isBedops);
	bool endMemoryProfile();
	bool calcMemoryStats();

};


#endif /* REGRESSTEST_H_ */
