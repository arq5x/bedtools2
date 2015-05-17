/*
 * sampleFile.h
 *
 *  Created on: Apr 29, 2015
 *      Author: nek3d
 */

#ifndef SAMPLEFILE_H_
#define SAMPLEFILE_H_

#include "ToolBase.h"
#include "ContextSample.h"

class SampleFile : public ToolBase {

public:
    SampleFile(ContextSample *context);
    virtual ~SampleFile();
	virtual bool init();
	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits) {}
	virtual void cleanupHits(RecordKeyVector &hits) {}
	virtual bool finalizeCalculations() { return true; }
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr);


protected:
	FileRecordMgr *_inputFile;
	vector<Record *> _samples;
	size_t _numSamples; //the number of samples we ultimately want
	size_t _numCurrSamples; //the number of samples kept so far.
	size_t _currRecordNum;
	int _seed;

	static const int DEFAULT_NUM_SAMPLES = 1000000;
	bool keepRecord(Record *record);
	bool strandComplies(const Record * record);

	virtual ContextSample *upCast(ContextBase *context) { return static_cast<ContextSample *>(context); }


};




#endif /* SAMPLEFILE_H_ */
