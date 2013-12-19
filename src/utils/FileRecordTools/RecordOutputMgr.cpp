/*
 * RecordOutputMgr.cpp
 *
 *  Created on: May 28, 2013
 *      Author: nek3d
 */

#include "RecordOutputMgr.h"
#include "Context.h"
#include "BlockMgr.h"
#include "Bed3Interval.h"
#include "Bed4Interval.h"
#include "BedGraphInterval.h"
#include "Bed5Interval.h"
#include "Bed6Interval.h"
#include "BedPlusInterval.h"
#include "Bed12Interval.h"
#include "BamRecord.h"
#include "VcfRecord.h"
#include "GffRecord.h"



#include <cstdio>


RecordOutputMgr::RecordOutputMgr()
: _context(NULL),
  _printable(true),
  _bamWriter(NULL),
  _currBlockList(NULL),
  _blockMgr(NULL)
{
	_blockMgr = new BlockMgr();
}

RecordOutputMgr::~RecordOutputMgr()
{
	if (_outBuf.size() > 0) {
		flush();
	}
	if (_bamWriter != NULL) {
		_bamWriter->Close();
		delete _bamWriter;
		_bamWriter = NULL;
	}
	delete _blockMgr;
	_blockMgr = NULL;

}

bool RecordOutputMgr::init(Context *context) {
	_context = context;
	_blockMgr->setContext(_context);
	if (_context->getOutputFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		//set-up BAM writer.
		_bamWriter = new BamTools::BamWriter();
		_bamWriter->SetCompressionMode(_context->getUncompressedBam() ?  BamTools::BamWriter::Uncompressed : BamTools::BamWriter::Compressed);

		_bamWriter->Open("stdout", _context->getHeader(_context->getBamHeaderAndRefIdx()).c_str(), _context->getReferences(_context->getBamHeaderAndRefIdx()));
	} else {
		//for everything but BAM, we'll copy output to an output buffer before printing.
		_outBuf.reserve(MAX_OUTBUF_SIZE);
	}
	if (_context->getProgram() == Context::INTERSECT) {
		if (_context->getAnyHit() || _context->getNoHit() || _context->getWriteCount()) {
			_printable = false;
		}
	}
	return true;
}

//void RecordOutputMgr::printHeader(const string &header)
//{
//	_outBuf.append(header);
//}

bool RecordOutputMgr::printKeyAndTerminate(RecordKeyList &keyList) {
	printBamType bamCode = printBamRecord(keyList);
	if (bamCode == BAM_AS_BAM) {
		return true;
	} else if (bamCode == NOT_BAM) {
		keyList.getKey()->print(_outBuf);
		return false;

	}
	//otherwise, it was BAM_AS_BED, and the key was printed.
	return false;

}

RecordOutputMgr::printBamType RecordOutputMgr::printBamRecord(RecordKeyList &keyList, bool bamOutputOnly)
{
	const Record *record = keyList.getKey();
	if (record->getType() == FileRecordTypeChecker::BAM_RECORD_TYPE) {
		if (_context->getOutputFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
			_bamWriter->SaveAlignment(static_cast<const BamRecord *>(record)->getAlignment());
			return BAM_AS_BAM;
		} else {
			if (!bamOutputOnly) {
				if (record->isUnmapped()) {
					record->printUnmapped(_outBuf);
				} else {
					static_cast<const BamRecord *>(record)->print(_outBuf, _currBlockList);
				}
			}
			return BAM_AS_BED;
		}
	}
	return NOT_BAM;
}

void RecordOutputMgr::printRecord(const Record *record)
{
	RecordKeyList keyList(record);
	printRecord(keyList);
}

void RecordOutputMgr::printRecord(RecordKeyList &keyList) {
	if (keyList.getKey()->getType() == FileRecordTypeChecker::BAM_RECORD_TYPE) {
		RecordKeyList blockList(keyList.getKey());
		bool deleteBlocks = false;
		_blockMgr->getBlocks(blockList, deleteBlocks);
		printRecord(keyList, &blockList);
		if (deleteBlocks) {
			_blockMgr->deleteBlocks(blockList);
		}
		return;
	}
    printRecord(keyList, NULL);

}

void RecordOutputMgr::printRecord(RecordKeyList &keyList, RecordKeyList *blockList)
{
	if (needsFlush()) {
		flush();
	}
	//The first time we print a record is when we print any header, because the header
	//hasn't been read from the query file until after the first record has also been read.
	if (_context->getPrintHeader()) {
		checkForHeader();
	}
	const_cast<Record *>(keyList.getKey())->undoZeroLength();
	_currBlockList = blockList;

	if (_context->getProgram() == Context::INTERSECT) {
		if (_printable) {
			if (keyList.empty()) {
				if (_context->getWriteAllOverlap()) {
					// -wao the user wants to force the reporting of 0 overlap
					if (printKeyAndTerminate(keyList)) {
						_currBlockList = NULL;
						return;
					}
					tab();
					null(false, true);
					tab();
					_outBuf.append('0');
					newline();
					if (needsFlush()) flush();
				}
				else if (_context->getLeftJoin()) {
					if (printKeyAndTerminate(keyList)) {
						_currBlockList = NULL;
						return;
					}
					tab();
					null(false, true);
					newline();
					if (needsFlush()) flush();
					_currBlockList = NULL;
					return;
				}
			} else {
				if (printBamRecord(keyList, true) == BAM_AS_BAM) {
					_currBlockList = NULL;
					return;
				}
				for (RecordKeyList::const_iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
					reportOverlapDetail(keyList.getKey(), iter->value());
				}
			}
		} else { // not printable
			reportOverlapSummary(keyList);
		}
		_currBlockList = NULL;
	} else if (_context->getProgram() == Context::SAMPLE) {
		if (!printKeyAndTerminate(keyList)) {
			newline();
		}
		_currBlockList = NULL;
		return;
	}
}

void RecordOutputMgr::checkForHeader() {
	if (_context->getProgram() == Context::INTERSECT) {
		if (_context->getPrintHeader()) {
			_outBuf.append(_context->getHeader(_context->getQueryFileIdx()));
		}
	} else if (_context->getProgram() == Context::SAMPLE) {
		if (_context->getPrintHeader()) {
			_outBuf.append(_context->getHeader(_context->getInputFileIdx()));
		}
	}
	_context->setPrintHeader(false);
	if (needsFlush()) flush();
}

void RecordOutputMgr::reportOverlapDetail(const Record *keyRecord, const Record *hitRecord)
{
	//get the max start and min end as strings.
	const_cast<Record *>(hitRecord)->undoZeroLength();


	const QuickString *startStr = NULL;
	const QuickString *endStr = NULL;
	int maxStart = 0;
	int minEnd = 0;

	int keyStart = keyRecord->getStartPos();
	int keyEnd = keyRecord->getEndPos();
	int hitStart = hitRecord->getStartPos();
	int hitEnd = hitRecord->getEndPos();

	if (  keyStart>= hitStart) {
		//the key start is after the hit start, but we need to check and make sure the hit end is at least after the keyStart.
		//The reason for this is that, in some rare cases, such as both the key and hit having been zero length intervals,
		//the normal process for intersection that allows us to simply report the maxStart and minEnd do not necessarily apply.
		if (hitEnd >= keyStart) {
			//this is ok. We have a normal intersection where the key comes after the hit.

			maxStart = keyStart;
			startStr = &(keyRecord->getStartPosStr());

			minEnd = min(keyEnd, hitEnd);
			endStr = keyRecord->getEndPos() < hitRecord->getEndPos() ? &(keyRecord->getEndPosStr()) : &(hitRecord->getEndPosStr());

		} else {
			//this is the weird case of not a "real" intersection. The keyStart is greater than the hitEnd. So just report the key as is.
			maxStart = keyStart;
			minEnd = keyEnd;
			startStr = &(keyRecord->getStartPosStr());
			endStr = &(keyRecord->getEndPosStr());
		}

	} else {
		//all of the above, but backwards. keyStart is before hitStart.
		if (keyEnd >= hitStart) {
			//normal intersection, key first
			maxStart = hitStart;
			startStr = &(hitRecord->getStartPosStr());
			minEnd = min(keyEnd, hitEnd);
			endStr = keyRecord->getEndPos() < hitRecord->getEndPos() ? &(keyRecord->getEndPosStr()) : &(hitRecord->getEndPosStr());
		} else {
			//this is the weird case of not a "real" intersection. The hitStart is greater than the keyEnd. So just report the hit as is.
			maxStart = hitStart;
			minEnd = hitEnd;
			startStr = &(hitRecord->getStartPosStr());
			endStr = &(hitRecord->getEndPosStr());

		}
	}

//	const QuickString &startStr = keyRecord->getStartPos() > hitRecord->getStartPos() ? keyRecord->getStartPosStr() : hitRecord->getStartPosStr();
//	const QuickString &endStr = keyRecord->getEndPos() < hitRecord->getEndPos() ? keyRecord->getEndPosStr() : hitRecord->getEndPosStr();
//
//	int maxStart = max(keyRecord->getStartPos(), hitRecord->getStartPos());
//	int minEnd = min(keyRecord->getEndPos(), hitRecord->getEndPos());
//

	if (!_context->getWriteA() && !_context->getWriteB() && !_context->getWriteOverlap() && !_context->getLeftJoin()) {
		printKey(keyRecord, *startStr, *endStr);
		newline();
		if (needsFlush()) flush();
	}
	else if ((_context->getWriteA() && _context->getWriteB()) || _context->getLeftJoin()) {
		printKey(keyRecord);
		tab();
		hitRecord->print(_outBuf);
		newline();
		if (needsFlush()) flush();
	}
	else if (_context->getWriteA()) {
		printKey(keyRecord);
		newline();
		if (needsFlush()) flush();
	}
	else if (_context->getWriteB()) {
		printKey(keyRecord, *startStr, *endStr);
		tab();
		hitRecord->print(_outBuf);
		newline();
		if (needsFlush()) flush();
	}
	else if (_context->getWriteOverlap()) {
		int printOverlapBases = max(0, minEnd-maxStart);
		printKey(keyRecord);
		tab();
		hitRecord->print(_outBuf);
		tab();
		int2str(printOverlapBases, _outBuf, true);
		newline();
		if (needsFlush()) flush();
	}
}

void RecordOutputMgr::reportOverlapSummary(RecordKeyList &keyList)
{
	int numOverlapsFound = (int)keyList.size();
	if (_context->getAnyHit() && numOverlapsFound > 0) {
		if (printKeyAndTerminate(keyList)) {
			return;
		}
		newline();
		if (needsFlush()) flush();
	} else if (_context->getWriteCount()) {
		if (printKeyAndTerminate(keyList)) {
			return;
		}
		tab();
		int2str(numOverlapsFound, _outBuf, true);
		newline();
		if (needsFlush()) flush();
	} else if (_context->getNoHit() && numOverlapsFound == 0) {
		if (printKeyAndTerminate(keyList)) {
			return;
		}
		newline();
		if (needsFlush()) flush();
	}
}


void RecordOutputMgr::null(bool queryType, bool dbType)
{
	FileRecordTypeChecker::RECORD_TYPE recordType = FileRecordTypeChecker::UNKNOWN_RECORD_TYPE;
	if (_context->getProgram() == Context::INTERSECT) {
		if (queryType) {
			recordType = _context->getQueryRecordType();
		} else if (dbType) {
			recordType = _context->getDatabaseRecordType();
		}
	} else if (_context->getProgram() == Context::SAMPLE) {
		recordType = _context->getInputRecordType();
	}
	//This is kind of a hack. Need an instance of the correct class of record in order to call it's printNull method.
	Record *dummyRecord = NULL;

	switch (recordType) {
	case FileRecordTypeChecker::BED3_RECORD_TYPE:
		dummyRecord = new Bed3Interval();
		break;
	case FileRecordTypeChecker::BED4_RECORD_TYPE:
		dummyRecord = new Bed4Interval();
		break;
	case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		dummyRecord = new BedGraphInterval();
		break;
	case FileRecordTypeChecker::BED5_RECORD_TYPE:
		dummyRecord = new Bed5Interval();
		break;
	case FileRecordTypeChecker::BED6_RECORD_TYPE:
		dummyRecord = new Bed6Interval();
		break;
	case FileRecordTypeChecker::BED12_RECORD_TYPE:
		dummyRecord = new Bed12Interval();
		break;
	case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		dummyRecord = new BedPlusInterval();
		(static_cast<BedPlusInterval *>(dummyRecord))->setNumPrintFields(_context->getMaxNumDatabaseFields());
		break;
	case FileRecordTypeChecker::VCF_RECORD_TYPE:
		dummyRecord = new VcfRecord();
		(static_cast<VcfRecord *>(dummyRecord))->setNumPrintFields(_context->getMaxNumDatabaseFields());
		break;
	case FileRecordTypeChecker::BAM_RECORD_TYPE:
		dummyRecord = new BamRecord();
		break;
	case FileRecordTypeChecker::GFF_RECORD_TYPE:
		dummyRecord = new GffRecord();
		(static_cast<GffRecord *>(dummyRecord))->setNumFields(_context->getMaxNumDatabaseFields());
		break;
	default:
		break;
	}

	dummyRecord->printNull(_outBuf);
	delete dummyRecord;

}

void RecordOutputMgr::printKey(const Record *key, const QuickString & start, const QuickString & end)
{
	if (key->getType() != FileRecordTypeChecker::BAM_RECORD_TYPE) {
		key->print(_outBuf, start, end);
	} else {
		static_cast<const BamRecord *>(key)->print(_outBuf, start, end, _currBlockList);
	}
}

void RecordOutputMgr::printKey(const Record *key)
{
	if (key->getType() != FileRecordTypeChecker::BAM_RECORD_TYPE) {
		key->print(_outBuf);
	} else {
		static_cast<const BamRecord *>(key)->print(_outBuf, _currBlockList);
	}
}

void RecordOutputMgr::flush() {
	fwrite(_outBuf.c_str(), 1, _outBuf.size(), stdout);
	_outBuf.clear();
}
