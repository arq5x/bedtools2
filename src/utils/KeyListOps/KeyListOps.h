/*
 * KeyListOps.h
 *
 *  Created on: Feb 24, 2014
 *      Author: nek3d
 */

#ifndef KEYLISTOPS_H_
#define KEYLISTOPS_H_

#include "KeyListOpsMethods.h"
#include "FileRecordTypeChecker.h"

class FileRecordMgr;

//print help message
void KeyListOpsHelp();

class KeyListOps {
public:

	KeyListOps();

	void setColumns(const QuickString &columns) { _columns = columns; }
	void addColumns(const QuickString &newCols) {
		if (!_columns.empty()) _columns += ",";
		_columns += newCols;
	}
	void setOperations(const QuickString & operation) { _operations = operation; }
	void addOperations(const QuickString &newOps) {
		if (!_operations.empty()) _operations += ",";
		_operations += newOps;
	}

	void setNullValue(const QuickString & nullValue) { _methods.setNullValue(nullValue); }
	void setDelimStr(const QuickString & delimStr) { _methods.setDelimStr(delimStr); }

	const QuickString &getColumns() { return _columns; }
	const QuickString &getOperations() { return _operations; }
	const QuickString &getNullValue() { return _methods.getNullValue(); }
	const QuickString &getDelimStr() { return _methods.getDelimStr(); }

	void setKeyList(RecordKeyVector *keyList) { _methods.setKeyList(keyList); }

	typedef enum { SUM, MEAN, STDDEV, SAMPLE_STDDEV, MEDIAN, MODE, ANTIMODE, MIN, MAX, ABSMIN, ABSMAX, COUNT, DISTINCT, COUNT_DISTINCT,
    	DISTINCT_ONLY, DISTINCT_SORT_NUM, DISTINCT_SORT_NUM_DESC, COLLAPSE, CONCAT, FREQ_ASC, FREQ_DESC, FIRST, LAST, INVALID } OP_TYPES;

	void setDBfileType(FileRecordTypeChecker::FILE_TYPE type) { _dbFileType = type; }
	bool isValidColumnOps(FileRecordMgr *dbFile);

	const QuickString &getOpVals(RecordKeyVector &hits);
	void setPrecision(int val) { _precision = val; }

private:
    void init();
    FileRecordTypeChecker::FILE_TYPE _dbFileType;

    QuickString _operations;
    QuickString _columns;

	KeyListOpsMethods _methods;
	map<QuickString, OP_TYPES> _opCodes;
	map<OP_TYPES, bool> _isNumericOp;

    typedef vector<pair<int, OP_TYPES> > colOpsType;
    colOpsType _colOps;
    QuickString _outVals;

    QuickString _formatStr;
    int _precision;

    static const int DEFAULT_PRECISION = 10;
    OP_TYPES getOpCode(const QuickString &operation) const;
    bool isNumericOp(OP_TYPES op) const;
    bool isNumericOp(const QuickString &op) const;
    const QuickString &format(double val);

};

#endif /* KEYLISTOPS_H_ */
