/*
 * SampleFile.h
 *
 *  Created on: Nov 18, 2013
 *      Author: nek3d
 */

#ifndef SAMPLEFILE_H_
#define SAMPLEFILE_H_

using namespace std;

#include "ContextSample.h"
#include "Record.h"
#include <vector>

class FileRecordMgr;
class Context;
class RecordOutputMgr;

class SampleFile {
public:
	SampleFile(ContextSample *context);
	~SampleFile();
	bool takeSample();

private:
	ContextSample *_context;
	FileRecordMgr *_inputFile;
	RecordOutputMgr *_outputMgr;
	vector<Record *> _samples;
	size_t _numSamples; //the number of samples we ultimately want
	size_t _numCurrSamples; //the number of samples kept so far.
	size_t _currRecordNum;
	int _seed;

	static const int DEFAULT_NUM_SAMPLES = 1000000;
	bool keepRecord(Record *record);
	bool strandComplies(const Record * record);

};

#endif /* SAMPLEFILE_H_ */
