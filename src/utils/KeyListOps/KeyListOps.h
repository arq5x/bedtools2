/*
 * KeyListOps.h
 *
 *  Created on: Feb 24, 2014
 *      Author: nek3d
 */

#ifndef KEYLISTOPS_H_
#define KEYLISTOPS_H_

#include "KeyListOpsMethods.h"

class FileRecordMgr;

class KeyListOps {
public:

	KeyListOps();

	void setColumns(const QuickString &columns) { _columns = columns; }
	void setOperations(const QuickString & operation) { _operations = operation; }
	void setNullValue(const QuickString & nullValue) { _methods.setNullValue(nullValue); }
	void setDelimStr(const QuickString & delimStr) { _methods.setDelimStr(delimStr); }

	void setKeyList(RecordKeyList *keyList) { _methods.setKeyList(keyList); }

	typedef enum { SUM, MEAN, STDDEV, SAMPLE_STDDEV, MEDIAN, MODE, ANTIMODE, MIN, MAX, ABSMIN, ABSMAX, COUNT, DISTINCT, COUNT_DISTINCT,
    	DISTINCT_ONLY, COLLAPSE, CONCAT, FREQ_ASC, FREQ_DESC, FIRST, LAST, INVALID } OP_TYPES;

	bool isValidColumnOps(FileRecordMgr *dbFile);

	const QuickString &getOpVals(RecordKeyList &hits);

private:
    void init();

    QuickString _operations;
    QuickString _columns;

	KeyListOpsMethods _methods;
	map<QuickString, OP_TYPES> _opCodes;
	map<OP_TYPES, bool> _isNumericOp;

    typedef vector<pair<int, OP_TYPES> > colOpsType;
    colOpsType _colOps;
    QuickString _outVals;

    OP_TYPES getOpCode(const QuickString &operation) const;
    bool isNumericOp(OP_TYPES op) const;
    bool isNumericOp(const QuickString &op) const;

};

#endif /* KEYLISTOPS_H_ */
