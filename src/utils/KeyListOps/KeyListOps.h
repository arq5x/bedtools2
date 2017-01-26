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

	void setColumns(const string &columns) { _columns = columns; }
	void addColumns(const string &newCols) {
		if (!_columns.empty()) _columns += ",";
		_columns += newCols;
	}
	void setOperations(const string & operation) { _operations = operation; }
	void addOperations(const string &newOps) {
		if (!_operations.empty()) _operations += ",";
		_operations += newOps;
	}

	void setNullValue(const string & nullValue) { _methods.setNullValue(nullValue); }
	void setDelimStr(const string & delimStr) { _methods.setDelimStr(delimStr); }

	const string &getColumns() { return _columns; }
	const string &getOperations() { return _operations; }
	const string &getNullValue() { return _methods.getNullValue(); }
	const string &getDelimStr() { return _methods.getDelimStr(); }

	void setKeyList(RecordKeyVector *keyList) { _methods.setKeyList(keyList); }

	typedef enum { SUM, MEAN, STDDEV, SAMPLE_STDDEV, MEDIAN, MODE, ANTIMODE, MIN, MAX, ABSMIN, ABSMAX, COUNT, DISTINCT, COUNT_DISTINCT,
    	DISTINCT_ONLY, DISTINCT_SORT_NUM, DISTINCT_SORT_NUM_DESC, COLLAPSE, CONCAT, FREQ_ASC, FREQ_DESC, FIRST, LAST, INVALID } OP_TYPES;

	void setDBfileType(FileRecordTypeChecker::FILE_TYPE type) { _dbFileType = type; }
	bool isValidColumnOps(FileRecordMgr *dbFile);

	const string &getOpVals(RecordKeyVector &hits);
	void setPrecision(int val) { _precision = val; }

private:
    void init();
    FileRecordTypeChecker::FILE_TYPE _dbFileType;

    string _operations;
    string _columns;

	KeyListOpsMethods _methods;
	map<string, OP_TYPES> _opCodes;
	map<OP_TYPES, bool> _isNumericOp;

    typedef vector<pair<int, OP_TYPES> > colOpsType;
    colOpsType _colOps;
    string _outVals;

    string _formatStr;
    int _precision;

    static const int DEFAULT_PRECISION = 10;
    OP_TYPES getOpCode(const string &operation) const;
    bool isNumericOp(OP_TYPES op) const;
    bool isNumericOp(const string &op) const;
    const string &format(double val);

};

#endif /* KEYLISTOPS_H_ */
