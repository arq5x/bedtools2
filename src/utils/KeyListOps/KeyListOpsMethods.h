/*
 * KeyListOpsMethods.h
 *
 *  Created on: Feb 6, 2014
 *      Author: nek3d
 */

#ifndef KEYLISTOPSMETHODS_H_
#define KEYLISTOPSMETHODS_H_

#include <map>
#include <utility> //for pair
#include "string.h"
#include <stdint.h>
#include "RecordKeyVector.h"

using namespace std;

class KeyListOpsMethods {
public:
	KeyListOpsMethods();
	KeyListOpsMethods(RecordKeyVector *keyList, int column = 1);
	~KeyListOpsMethods();

	void setIsBam(bool isBam) { _isBam = isBam; }
	void setKeyList(RecordKeyVector *keyList) { _keyList = keyList; }
	void setColumn(int col) { _column = col; }
	void setNullValue(const string & nullVal) { _nullVal = nullVal; }
	const string &getNullValue() const { return _nullVal; }
	void setDelimStr(const string &delimStr) { _delimStr = delimStr; }
	const string &getDelimStr() const { return _delimStr; }

    // return the total of the values in the vector
    double getSum();
    // return the average value in the vector
    double getMean();
     // return the standard deviation
    double getStddev();
    // return the sample standard deviation
    double getSampleStddev();
    // return the median value in the vector
    double getMedian();
    // return the most common value in the vector
    const string &getMode();
    // return the least common value in the vector
    const string &getAntiMode();
    // return the minimum element of the vector
    double getMin();
    // return the maximum element of the vector
    double getMax();
    // return the minimum absolute value of the vector
    double getAbsMin();
    // return the maximum absolute value of the vector
    double getAbsMax();
    // return the count of element in the vector
    uint32_t getCount();
    // return a the count of _unique_ elements in the vector
    uint32_t getCountDistinct();
    // return only those elements that occur once
    const string &getDistinctOnly();
    // as distinct, but sorted numerically.
    const string &getDistinctSortNum(bool ascending = true);
    // as distinct, but sorted numerically, descending


     // return a delimiter-separated list of elements
    const string & getCollapse(const string & delimiter = ",");
    // return a concatenation of all elements in the vector
    const string & getConcat();
    // return a comma-separated list of the _unique_ elements
    const string & getDistinct();
    // return a histogram of values and their freqs. in desc. order of frequency
    const string & getFreqDesc();
    // return a histogram of values and their freqs. in asc. order of frequency
    const string & getFreqAsc();
    // return the first value in the list
    const string & getFirst();
    // return the last value in the list
    const string & getLast();

    bool nonNumErrFlagSet() const { return _nonNumErrFlag; }
    const string &getErrMsg() const { return _errMsg; }
    void resetNonNumErrFlag() {
    	_nonNumErrFlag = false;
    	_errMsg.clear();
    }

private:
	RecordKeyVector *_keyList;
	int _column;
	string _nullVal;
	string _delimStr;
	string _retStr;

	RecordKeyVector _nullKeyList; //this has to exist just so we can initialize _iter, below.
	RecordKeyVector::iterator_type _iter;

	// Some methods need to put values into a vector, mostly for sorting.
	vector<double> _numArray;
	vector<string> _qsArray;

	typedef map<string, int> freqMapType;
	freqMapType _freqMap;
	freqMapType::iterator _freqIter;

	typedef enum { UNSORTED, ASC, DESC} SORT_TYPE;

	bool _nonNumErrFlag;
	string _errMsg;
	bool _isBam;

	typedef multimap<int, string, less<int> > histAscType;
	typedef multimap<int, string, greater<int> > histDescType;
	void init();
	const string &getColVal();
	double getColValNum();
	bool empty() { return _keyList->empty(); }
	void begin() { _iter = _keyList->begin(); }
	bool end() { return _iter == _keyList->end(); }
	void next() { _iter = _keyList->next(); }
	void toArray(bool useNum, SORT_TYPE sortVal = UNSORTED);
	void sortArray(bool useNum, bool ascOrder);
	void makeFreqMap();


};


#endif /* KEYLISTOPSMETHODS_H_ */
