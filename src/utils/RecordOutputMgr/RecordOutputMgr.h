/*
 * RecordOutputMgr.h
 *
 *  Created on: May 28, 2013
 *      Author: nek3d
 */

#ifndef RECORDOUTPUTMGR_H_
#define RECORDOUTPUTMGR_H_

#include "ContextBase.h"
#include "RecordKeyVector.h"
#include "api/BamWriter.h"

using namespace std;

class BlockMgr;

class RecordOutputMgr {
public:
	RecordOutputMgr();
	~RecordOutputMgr();

	//The init method must be called after all the input files are open.
	void init(ContextBase *context);

	void printRecord(Record *record);
	void printRecord(RecordKeyVector &keyList);
	void printRecord(Record *record, const string & value);
	void checkForHeader();

	void printClosest(RecordKeyVector &keyList, const vector<CHRPOS> *dists = NULL);

	void tab() { _outBuf.append("\t"); }
	void newline() { _outBuf.append("\n"); }

private:
	typedef enum { NOT_BAM, BAM_AS_BAM, BAM_AS_BED} printBamType;

	ContextBase *_context;
	bool _printable;
	BamTools::BamWriter *_bamWriter;
	RecordKeyVector *_currBamBlockList;

	string _outBuf;

	BlockMgr *_bamBlockMgr;
	string _afterVal; //to store values to be printed after record, such as column operations.
	//some helper functions to neaten the code.
	void null(bool queryType, bool dbType);

	void printRecord(RecordKeyVector &keyList, RecordKeyVector *blockList);
	void printKey(const Record *key);
	void printKey(const Record *key, const string & start, const string & end);
	void printKey(const Record *key, CHRPOS start, CHRPOS end);
	void addDbFileId(int fileId);
	bool printKeyAndTerminate(RecordKeyVector &keyList);
	printBamType printBamRecord(RecordKeyVector &keyList, bool bamOutputOnly = false);
	void reportOverlapDetail(const Record *keyRecord, const Record *hitRecord, int hitIdx = 0);
	void reportOverlapSummary(RecordKeyVector &keyList);

	static const unsigned int MAX_OUTBUF_SIZE = 16384; //16 K

	// If we are using buffered output, only flush the output buffer if it's least
	// 90% full. If we're not using buffered output, flush if it's not empty
	bool needsFlush() const {
		return ((_context->getUseBufferedOutput() &&_outBuf.size() >= MAX_OUTBUF_SIZE *.9) ||
				(!_context->getUseBufferedOutput() && !_outBuf.empty()));
	}
	void flush();
};

#endif /* RECORDOUTPUTMGR_H_ */
