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


class RecordOutputMgr {
public:
	RecordOutputMgr();
	~RecordOutputMgr();

	//The init method must be called after all the input files are open.
	bool init(Context *context);

	typedef enum { NOT_BAM, BAM_AS_BAM, BAM_AS_BED} printBamType;

	void printHeader(const string &);
	bool printKeyAndTerminate(RecordKeyList &keyList);
	printBamType printBamRecord(RecordKeyList &keyList, bool bamOutputOnly = false);
	void printRecord(RecordKeyList &keyList, RecordKeyList *blockList = NULL);
	void reportOverlapDetail(const Record *keyRecord, const Record *hitRecord);
	void reportOverlapSummary(RecordKeyList &keyList);

private:
	Context *_context;
	bool _printable;
	BamTools::BamWriter *_bamWriter;
	RecordKeyList *_currBlockList;

	QuickString _outBuf;

	//some helper functions to neaten the code.
	void tab() { _outBuf.append('\t'); }
	void newline() { _outBuf.append('\n'); }
	void null(bool queryType, bool dbType);

	void printKey(const Record *key);
	void printKey(const Record *key, const QuickString & start, const QuickString & end);
	static const unsigned int MAX_OUTBUF_SIZE = 33554432; //32 Mb
	bool needsFlush() const { return _outBuf.size() >= MAX_OUTBUF_SIZE *.9; }
	void flush();
};

#endif /* RECORDOUTPUTMGR_H_ */
