/*
 * KeyListOps.cpp
 *
 *  Created on: Feb 24, 2014
 *      Author: nek3d
 */
#include "KeyListOps.h"
#include "FileRecordMgr.h"
#include <cmath> //for std::isnan
#include <sstream>
#include <iomanip>

KeyListOps::KeyListOps():
_dbFileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE)
{
	_opCodes["sum"] = SUM;
	_opCodes["mean"] = MEAN;
	_opCodes["stdev"] = STDDEV;
	_opCodes["sstdev"] = SAMPLE_STDDEV;
	_opCodes["median"] = MEDIAN;
	_opCodes["mode"] = MODE;
	_opCodes["antimode"] = ANTIMODE;
	_opCodes["min"] = MIN;
	_opCodes["max"] = MAX;
	_opCodes["absmin"] = ABSMIN;
	_opCodes["absmax"] = ABSMAX;
	_opCodes["count"] = COUNT;
	_opCodes["distinct"] = DISTINCT;
	_opCodes["count_distinct"] = COUNT_DISTINCT;
	_opCodes["distinct_only"] = DISTINCT_ONLY;
	_opCodes["distinct_sort_num"] = DISTINCT_SORT_NUM;
	_opCodes["distinct_sort_num_desc"] = DISTINCT_SORT_NUM_DESC;

	_opCodes["collapse"] = COLLAPSE;
	_opCodes["concat"] = CONCAT;
	_opCodes["freqasc"] = FREQ_ASC;
	_opCodes["freqdesc"] = FREQ_DESC;
	_opCodes["first"] = FIRST;
	_opCodes["last"] = LAST;

	_isNumericOp[SUM] = true;
	_isNumericOp[MEAN] = true;
	_isNumericOp[STDDEV] = true;
	_isNumericOp[MEDIAN] = true;
	_isNumericOp[MODE] = false;
	_isNumericOp[ANTIMODE] = false;
	_isNumericOp[MIN] = true;
	_isNumericOp[MAX] = true;
	_isNumericOp[ABSMIN] = true;
	_isNumericOp[COUNT] = false;
	_isNumericOp[DISTINCT] = false;
	_isNumericOp[COUNT_DISTINCT] = false;
	_isNumericOp[DISTINCT_ONLY] = false;
	_isNumericOp[COLLAPSE] = false;
	_isNumericOp[CONCAT] = false;
	_isNumericOp[FREQ_ASC] = false;
	_isNumericOp[FREQ_DESC] = false;
	_isNumericOp[FIRST] = false;
	_isNumericOp[LAST] = false;

	_methods.setDelimStr(",");
	_methods.setNullValue(".");

	// default to BED score column
	_columns = "5";
	// default to "sum"
	_operations = "sum";
	_precision = DEFAULT_PRECISION;

}

bool KeyListOps::isNumericOp(OP_TYPES op) const {
	map<OP_TYPES, bool>::const_iterator iter = _isNumericOp.find(op);
	return (iter == _isNumericOp.end() ? false : iter->second);
}

bool KeyListOps::isNumericOp(const string &op) const {
	return isNumericOp(getOpCode(op));
}

KeyListOps::OP_TYPES KeyListOps::getOpCode(const string &operation) const {
	//If the operation does not exist, return INVALID.
	//otherwise, return code for given operation.
	map<string, OP_TYPES>::const_iterator iter = _opCodes.find(operation);
	if (iter == _opCodes.end()) {
		return INVALID;
	}
	return iter->second;
}


bool KeyListOps::isValidColumnOps(FileRecordMgr *dbFile) {

	//get the strings from context containing the comma-delimited lists of columns
	//and operations. Split both of these into vectors. Get the operation code
	//for each operation string. Finally, make a vector of pairs, where the first
	//member of each pair is a column number, and the second member is the code for the
	//operation to perform on that column.

    Tokenizer colTokens;
    Tokenizer opsTokens;

    int numCols = colTokens.tokenize(_columns, ',');
	int numOps = opsTokens.tokenize(_operations, ',');

	if (numOps < 1 || numCols < 1) {
		 cerr << endl << "*****" << endl
		             << "***** ERROR: There must be at least one column and at least one operation named." << endl;
		 return false;
	}
	if (numOps > 1 && numCols > 1 && numCols != numOps) {
		 cerr << endl << "*****" << endl
		             << "***** ERROR: There are " << numCols <<" columns given, but there are " << numOps << " operations." << endl;
		cerr << "\tPlease provide either a single operation that will be applied to all listed columns, " << endl;
		cerr << "\ta single column to which all operations will be applied," << endl;
		cerr << "\tor an operation for each column." << endl;
		return false;
	}
	int loop = max(numCols, numOps);

	// If there is only one column, all ops are performed on it.
	// Otherwise, if there is only op, it is performed on all columns.
	// Besides that, ops are performed on columns in their respective
	// ordering.

	for (int i=0; i < loop; i++) {
		int col = (int)str2chrPos(colTokens.getElem(numCols > 1 ? i : 0));

		//check that the column number is valid
		if (col < 1 || col > dbFile->getNumFields()) {
			 cerr << endl << "*****" << endl  << "***** ERROR: Requested column " << col << ", but database file "
					 << dbFile->getFileName() << " only has fields 1 - " << dbFile->getNumFields() << "." << endl;
			 return false;
		}
		const string &operation = opsTokens.getElem(numOps > 1 ? i : 0);
		OP_TYPES opCode = getOpCode(operation);
		if (opCode == INVALID) {
			cerr << endl << "*****" << endl
								 << "***** ERROR: " << operation << " is not a valid operation. " << endl;
			return false;
		}
		_colOps.push_back(pair<int, OP_TYPES>(col, opCode));
	}

	//lastly, if the file is BAM, and they asked for column 2, which is the
	//flags field, then for now we have to throw an error, as the flag field
	//is currently not supported.
	if (_dbFileType == FileRecordTypeChecker::BAM_FILE_TYPE) {
		//also, tell the methods class we're dealing with BAM.
		_methods.setIsBam(true);
		for (size_t i = 0; i < _colOps.size(); i++) {
			if (_colOps[i].first == 2) {
				cerr << endl << "*****" << endl << "***** ERROR: Requested column 2 of a BAM file, which is the Flags field." << endl;
				cerr << "             We currently do not support this, but may in future versions." << endl;
				return false;
			}
		}
	}

    return true;
}

const string & KeyListOps::getOpVals(RecordKeyVector &hits)
{
	//loop through all requested columns, and for each one, call the method needed
	//for the operation specified.
	_methods.setKeyList(&hits);
	_outVals.clear();
	double val = 0.0;
	for (int i=0; i < (int)_colOps.size(); i++) {
		ostringstream s;
		int col = _colOps[i].first;
		OP_TYPES opCode = _colOps[i].second;

		_methods.setColumn(col);
		switch (opCode) {
		case SUM:
			val = _methods.getSum();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case MEAN:
			val = _methods.getMean();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case STDDEV:
			val = _methods.getStddev();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case SAMPLE_STDDEV:
			val = _methods.getSampleStddev();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case MEDIAN:
			val = _methods.getMedian();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case MODE:
			s << _methods.getMode();
			break;

		case ANTIMODE:
			s << _methods.getAntiMode();
			break;

		case MIN:
			val = _methods.getMin();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case MAX:
			val = _methods.getMax();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case ABSMIN:
			val = _methods.getAbsMin();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case ABSMAX:
			val = _methods.getAbsMax();
			if (std::isnan(val)) {
				s << _methods.getNullValue();
			} else {
				s << format(val);
			}
			break;

		case COUNT:
			s << _methods.getCount();
			break;

		case DISTINCT:
			s << _methods.getDistinct();
			break;

		case DISTINCT_SORT_NUM:
			s << _methods.getDistinctSortNum();
			break;

		case DISTINCT_SORT_NUM_DESC:
			s << _methods.getDistinctSortNum(false);
			break;

		case COUNT_DISTINCT:
			s << _methods.getCountDistinct();
			break;

		case DISTINCT_ONLY:
			s << _methods.getDistinctOnly();
			break;

		case COLLAPSE:
			s << _methods.getCollapse();
			break;

		case CONCAT:
			s << _methods.getConcat();
			break;

		case FREQ_ASC:
			s << _methods.getFreqAsc();
			break;

		case FREQ_DESC:
			s << _methods.getFreqDesc();
			break;

		case FIRST:
			s << _methods.getFirst();
			break;

		case LAST:
			s << _methods.getLast();
			break;

		case INVALID:
		default:
			// Any unrecognized operation should have been handled already in the context validation.
			// It's thus unnecessary to handle it here, but throw an error to help us know if future
			// refactoring or code changes accidentally bypass the validation phase.
			cerr << "ERROR: Invalid operation given for column " << col << ". Exiting..." << endl;
			break;
		}
		_outVals.append(s.str());
		//if this isn't the last column, add a tab.
		if (i < (int)_colOps.size() -1) {
			_outVals.append("\t");
		}
	}
	if (_methods.nonNumErrFlagSet()) {
		//asked for a numeric op on a column in which a non numeric value was found.
		cerr << _methods.getErrMsg() << endl;
		_methods.resetNonNumErrFlag();
	}
	return _outVals;
}

const string &KeyListOps::format(double val)
{
   std::stringstream strmBuf;
   strmBuf << std::setprecision (_precision) << val;
   _formatStr = strmBuf.str();
   return _formatStr;
}

void KeyListOpsHelp() {

    cerr << "\t-c\t"             << "Specify columns from the B file to map onto intervals in A." << endl;
    cerr                         << "\t\tDefault: 5." << endl;
    cerr						<<  "\t\tMultiple columns can be specified in a comma-delimited list." << endl << endl;

    cerr << "\t-o\t"             << "Specify the operation that should be applied to -c." << endl;
    cerr                         << "\t\tValid operations:" << endl;
    cerr                         << "\t\t    sum, min, max, absmin, absmax," << endl;
    cerr                         << "\t\t    mean, median, mode, antimode" << endl;
    cerr                         << "\t\t    stdev, sstdev" << endl;
    cerr                         << "\t\t    collapse (i.e., print a delimited list (duplicates allowed)), " << endl;
    cerr                         << "\t\t    distinct (i.e., print a delimited list (NO duplicates allowed)), " << endl;
    cerr                         << "\t\t    distinct_sort_num (as distinct, sorted numerically, ascending)," << endl;
    cerr                         << "\t\t    distinct_sort_num_desc (as distinct, sorted numerically, desscending)," << endl;
    cerr                         << "\t\t    distinct_only (delimited list of only unique values)," << endl;
    cerr                         << "\t\t    count" << endl;
    cerr                         << "\t\t    count_distinct (i.e., a count of the unique values in the column), " << endl;
    cerr                         << "\t\t    first (i.e., just the first value in the column), " << endl;
    cerr                         << "\t\t    last (i.e., just the last value in the column), " << endl;
    cerr                         << "\t\tDefault: sum" << endl;
    cerr						 << "\t\tMultiple operations can be specified in a comma-delimited list." << endl << endl;

    cerr						<< "\t\tIf there is only column, but multiple operations, all operations will be" << endl;
    cerr						<< "\t\tapplied on that column. Likewise, if there is only one operation, but" << endl;
    cerr						<< "\t\tmultiple columns, that operation will be applied to all columns." << endl;
    cerr						<< "\t\tOtherwise, the number of columns must match the the number of operations," << endl;
    cerr						<< "\t\tand will be applied in respective order." << endl;
    cerr						<< "\t\tE.g., \"-c 5,4,6 -o sum,mean,count\" will give the sum of column 5," << endl;
    cerr						<< "\t\tthe mean of column 4, and the count of column 6." << endl;
    cerr						<< "\t\tThe order of output columns will match the ordering given in the command." << endl << endl<<endl;

    cerr << "\t-delim\t"                 << "Specify a custom delimiter for the collapse operations." << endl;
    cerr                                 << "\t\t- Example: -delim \"|\"" << endl;
    cerr                                 << "\t\t- Default: \",\"." << endl << endl;

    cerr << "\t-prec\t"   		 << "Sets the decimal precision for output (Default: 5)" << endl << endl;
}
