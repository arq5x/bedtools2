/* groupBy.cpp
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#include "groupBy.h"
#include "Tokenizer.h"
#include "ParseTools.h"
#include "stringUtilities.h"
#include <utility>

GroupBy::GroupBy(ContextGroupBy *context)
: ToolBase(context),
  _queryFRM(NULL),
  _prevRecord(NULL)
{

}

GroupBy::~GroupBy()
{

}

bool GroupBy::init()
{
	Tokenizer groupColsTokens;
	groupColsTokens.tokenize(upCast(_context)->getGroupCols(), ',');
	int numElems = groupColsTokens.getNumTotalElems();
	for (int i=0; i < numElems; i++) {
		//if the item is a range, such as 3-5,
		//must split that as well.

		const string &elem = groupColsTokens.getElem(i);
		if (strchr(elem.c_str(), '-')) {
			Tokenizer rangeElems;
			rangeElems.tokenize(elem, '-');
			int startNum = (int)str2chrPos(rangeElems.getElem(0));
			int endNum = (int)str2chrPos(rangeElems.getElem(1));
			for (int i=startNum; i <= endNum; i++) {
				_groupCols.push_back(i);
			}
		} else {
			_groupCols.push_back(str2chrPos(elem));
		}
	}
	_queryFRM = _context->getFile(0);
	_prevFields.resize(_groupCols.size());

	_prevRecord = getNextRecord();
	return true;
}

bool GroupBy::findNext(RecordKeyVector &hits)
{
	//get one record.
	if (_prevRecord == NULL) {
		return false;
	}
	assignPrevFields();
	hits.setKey(_prevRecord);
	hits.push_back(_prevRecord); //key should also be part of group for calculations
	while (1) 
	{
		Record *newRecord = getNextRecord();
		if (newRecord == NULL) 
		{
			_prevRecord = NULL;
			break;
		} else if (canGroup(newRecord)) 
		{
			hits.push_back(newRecord);
		} 
		else 
		{
			_prevRecord = newRecord;
			break;
		}
	}
	return true;
}

void GroupBy::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{

	Record *rec = hits.getKey();
	const string &opVal  = _context->getColumnOpsVal(hits);
	if (upCast(_context)->printFullCols()) 
	{
		outputMgr->printRecord(rec, opVal);
	} 
	else 
	{
		string outBuf;
		for (int i = 0; i < (int)_groupCols.size(); i++) 
		{
			outBuf.append(rec->getField(_groupCols[i]));
			outBuf.append("\t");
		}
		outBuf.append(opVal);
		outputMgr->printRecord(NULL, outBuf);
	}
	cleanupHits(hits);
}

void GroupBy::cleanupHits(RecordKeyVector &hits)
{
	RecordKeyVector::iterator_type iter = hits.begin();
	for (; iter != hits.end(); iter = hits.next()) 
	{
		_queryFRM->deleteRecord(*iter);	
	}
	hits.clearAll();
}

Record *GroupBy::getNextRecord() {
	while (!_queryFRM->eof()) 
	{
		Record *queryRecord = _queryFRM->getNextRecord();
		if (queryRecord == NULL) 
		{
			continue;
		} 
		else 
		{
			return queryRecord;
		}
	}
	return NULL;
}

void GroupBy::assignPrevFields() {
	for (int i=0; i < (int)_prevFields.size(); i++) {
		_prevFields[i] = _prevRecord->getField(_groupCols[i]);
	}
}

bool GroupBy::canGroup(Record *newRecord) 
{
	for (int i = 0; i < (int)_groupCols.size(); i++) 
	{
		int fieldNum = _groupCols[i];
		const string &newField = newRecord->getField(fieldNum);
		const string &oldField = _prevFields[i];
		if (upCast(_context)->ignoreCase()) 
		{
			if (toLower(oldField) != toLower(newField)) return false;
		} 
		else 
		{
			if (oldField != newField) return false;
		}
	}
	return true;
}

