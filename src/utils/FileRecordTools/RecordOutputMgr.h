/*
 * RecordOutputMgr.h
 *
 *  Created on: May 28, 2013
 *      Author: nek3d
 */

#ifndef RECORDOUTPUTMGR_H_
#define RECORDOUTPUTMGR_H_

using namespace std;

#include "RecordKeyList.h"
#include "api/BamWriter.h"

class Context;
class BlockMgr;

class RecordOutputMgr {
public:
	RecordOutputMgr();
	~RecordOutputMgr();

	//The init method must be called after all the input files are open.
	bool init(Context *context);

//	void printHeader(const string &);
	void printRecord(const Record *record);
	void printRecord(RecordKeyList &keyList);

private:
	typedef enum { NOT_BAM, BAM_AS_BAM, BAM_AS_BED} printBamType;

	Context *_context;
	bool _printable;
	BamTools::BamWriter *_bamWriter;
	RecordKeyList *_currBlockList;

	QuickString _outBuf;
	BlockMgr *_blockMgr;
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
	void reportOverlapDetail(const Record *keyRecord, const Record *hitRecord);
	void reportOverlapSummary(RecordKeyList &keyList);

	static const unsigned int MAX_OUTBUF_SIZE = 16384; //16 K
	bool needsFlush() const { return _outBuf.size() >= MAX_OUTBUF_SIZE *.9; }
	void flush();
};

#endif /* RECORDOUTPUTMGR_H_ */
