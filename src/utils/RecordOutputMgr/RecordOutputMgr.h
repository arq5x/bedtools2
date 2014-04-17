/*
 * RecordOutputMgr.h
 *
 *  Created on: May 28, 2013
 *      Author: nek3d
 */

#ifndef RECORDOUTPUTMGR_H_
#define RECORDOUTPUTMGR_H_

using namespace std;

#include "ContextBase.h"
#include "RecordKeyList.h"
#include "api/BamWriter.h"

class BlockMgr;

class RecordOutputMgr {
public:
	RecordOutputMgr();
	~RecordOutputMgr();

	//The init method must be called after all the input files are open.
	void init(ContextBase *context);

	void printRecord(const Record *record);
	void printRecord(RecordKeyList &keyList);
	void printRecord(const Record *record, const QuickString & value);

	//where necessary, pass additional information about splits through the blockMgr.
	void setSplitInfo(const BlockMgr *blockMgr) { _splitInfo = blockMgr; }

private:
	typedef enum { NOT_BAM, BAM_AS_BAM, BAM_AS_BED} printBamType;

	ContextBase *_context;
	bool _printable;
	BamTools::BamWriter *_bamWriter;
	RecordKeyList *_currBamBlockList;

	QuickString _outBuf;

	//
	BlockMgr *_bamBlockMgr;
	const BlockMgr *_splitInfo;
	QuickString _afterVal; //to store values to be printed after record, such as column operations.
	//some helper functions to neaten the code.
	void tab() { _outBuf.append('\t'); }
	void newline() { _outBuf.append('\n'); }
	void null(bool queryType, bool dbType);

	void printRecord(RecordKeyList &keyList, RecordKeyList *blockList);
	void printKey(const Record *key);
	void printKey(const Record *key, const QuickString & start, const QuickString & end);
	bool printKeyAndTerminate(RecordKeyList &keyList);
	printBamType printBamRecord(RecordKeyList &keyList, bool bamOutputOnly = false);
	void checkForHeader();
	void reportOverlapDetail(const Record *keyRecord, const Record *hitRecord, int hitIdx = 0);
	void reportOverlapSummary(RecordKeyList &keyList);

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
