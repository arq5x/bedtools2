/*
 * KeyListOpsMethods.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: nek3d
 */

#include "KeyListOpsMethods.h"
#include <cmath>
#include <algorithm>
#include <limits.h>
#include "ParseTools.h" //to get the isNumeric function

KeyListOpsMethods::KeyListOpsMethods()
: _keyList(&_nullKeyList),
  _column(1),
  _nullVal("."),
  _delimStr(","),
  _iter(_nullKeyList.begin()),
  _nonNumErrFlag(false),
  _isBam(false)
{
}

KeyListOpsMethods::KeyListOpsMethods(RecordKeyVector *keyList, int column)
: _keyList(keyList),
  _column(column),
  _nullVal("."),
  _delimStr(","),
  _iter(keyList->begin())
{
}


KeyListOpsMethods::~KeyListOpsMethods() {

}

// return the total of the values in the vector
double KeyListOpsMethods::getSum() {
	if (empty()) return NAN;

	double theSum = 0.0;
	for (begin(); !end(); next()) {
		theSum += getColValNum();
	}
	return theSum;
}

// return the average value in the vector
double KeyListOpsMethods::getMean() {
	if (empty()) return NAN;

	return getSum() / (float)getCount();
}


 // return the standard deviation
double KeyListOpsMethods::getStddev() {
	if (empty()) return NAN;

	double avg = getMean();
	double squareDiffSum = 0.0;
	for (begin(); !end(); next()) {
		double val = getColValNum();
		double diff = val - avg;
		squareDiffSum += diff * diff;
	}
	return sqrt(squareDiffSum / (float)getCount());
}
// return the standard deviation
double KeyListOpsMethods::getSampleStddev() {
	if (empty()) return NAN;

	double avg = getMean();
	double squareDiffSum = 0.0;
	for (begin(); !end(); next()) {
		double val = getColValNum();
		double diff = val - avg;
		squareDiffSum += diff * diff;
	}
	return sqrt(squareDiffSum / ((float)getCount() - 1.0));
}

// return the median value in the vector
double KeyListOpsMethods::getMedian() {
	if (empty()) return NAN;

	//get sorted vector. if even number of elems, return middle val.
	//if odd, average of two.
	toArray(true, ASC);
	size_t count = getCount();
	if (count % 2) {
		//odd number of elements. Take middle one.
		return _numArray[count/2];
	} else {
		//even numnber of elements. Take average of middle 2.
		double sum = _numArray[count/2 -1] + _numArray[count/2];
		return sum / 2.0;
	}
}

// return the most common value in the vector
const string &KeyListOpsMethods::getMode() {
	if (empty()) return _nullVal;

	makeFreqMap();

	//now pass through the freq map and keep track of which key has the highest occurance.
	freqMapType::iterator maxIter = _freqMap.begin();
	int maxVal = 0;
	for (; _freqIter != _freqMap.end(); _freqIter++) {
		if (_freqIter->second > maxVal) {
			maxIter = _freqIter;
			maxVal = _freqIter->second;
		}
	}
	_retStr = maxIter->first;
	return _retStr;
}
// return the least common value in the vector
const string &KeyListOpsMethods::getAntiMode() {
	if (empty()) return _nullVal;

	makeFreqMap();

	//now pass through the freq map and keep track of which key has the highest occurance.
	freqMapType::iterator minIter = _freqMap.begin();
	int minVal = INT_MAX;
	for (; _freqIter != _freqMap.end(); _freqIter++) {
		if (_freqIter->second < minVal) {
			minIter = _freqIter;
			minVal = _freqIter->second;
		}
	}
	_retStr =  minIter->first;
	return _retStr;
}
// return the minimum element of the vector
double KeyListOpsMethods::getMin() {
	if (empty()) return NAN;

	begin();
	double minVal = getColValNum();
	for (; !end(); next()) {
		double currVal = getColValNum();
		minVal = (currVal < minVal) ? currVal : minVal;
	}
	return  minVal;
}

// return the maximum element of the vector
double KeyListOpsMethods::getMax() {
	if (empty()) return NAN;

	begin();
	double maxVal = getColValNum();
	for (; !end(); next()) {
		double currVal = getColValNum();
		maxVal = (currVal > maxVal) ? currVal : maxVal;
	}
	return maxVal;
}

// return the minimum absolute value of the vector
double KeyListOpsMethods::getAbsMin() {
	if (empty()) return NAN;

	begin();
	double minVal = abs(getColValNum());
	for (; !end(); next()) {
		double currVal = abs(getColValNum());
		minVal = (currVal < minVal) ? currVal : minVal;
	}
	return minVal;
}
// return the maximum absolute value of the vector
double KeyListOpsMethods::getAbsMax() {
	if (empty()) return NAN;

	begin();
	double maxVal = abs(getColValNum());
	for (; !end(); next()) {
		double currVal = abs(getColValNum());
		maxVal = (currVal > maxVal) ? currVal : maxVal;
	}
	return maxVal;
}
// return the count of element in the vector
uint32_t KeyListOpsMethods::getCount() {
	return _keyList->size();
}
// return a delimited list of the unique elements
const string &KeyListOpsMethods::getDistinct() {
	if (empty()) return _nullVal;
	// separated list of unique values. If something repeats, only report once.
	makeFreqMap();
	_retStr.clear();
	for (; _freqIter != _freqMap.end(); _freqIter++) {
		if (_freqIter != _freqMap.begin()) _retStr += _delimStr;
		_retStr.append(_freqIter->first);
	}
	return _retStr;
}

const string &KeyListOpsMethods::getDistinctOnly() {
	if (empty()) return _nullVal;
	// separated list of unique values. If something repeats, don't report.
	makeFreqMap();
	_retStr.clear();
	for (; _freqIter != _freqMap.end(); _freqIter++) {
		if (_freqIter->second > 1) continue;
		if (_freqIter != _freqMap.begin()) _retStr += _delimStr;
		_retStr.append(_freqIter->first);
	}
	return _retStr;
}

const string &KeyListOpsMethods::getDistinctSortNum(bool asc) {
	if (empty()) return _nullVal;
	
	toArray(true, asc ? ASC : DESC);
	vector<double>::iterator endIter = std::unique(_numArray.begin(), _numArray.end());

	_retStr.clear();
	ostringstream s;
	for (vector<double>::iterator iter = _numArray.begin(); iter != endIter; iter++) {
		if (iter != _numArray.begin()) s << _delimStr;
		s << *iter;
	}
	_retStr.append(s.str());
	return  _retStr;

}


// return a the count of _unique_ elements in the vector
uint32_t KeyListOpsMethods::getCountDistinct() {
	if (empty()) return 0;

	makeFreqMap();
	return _freqMap.size();
}

// return a delimiter-separated list of elements
const string &KeyListOpsMethods::getCollapse(const string &delimiter) {
	if (empty()) return _nullVal;

	//just put all items in one big separated list.
	_retStr.clear();
	int i=0;
	for (begin(); !end(); next()) {
		if (i > 0) _retStr += _delimStr;
		_retStr.append(getColVal());
		i++;
	}
	return _retStr;
}

// return a concatenation of all elements in the vector
const string &KeyListOpsMethods::getConcat() {
	if (empty()) return _nullVal;

	//like collapse but w/o commas. Just a true concat of all vals.
	//just swap out the delimChar with '' and call collapse, then
	//restore the delimChar.
	string oldDelimStr(_delimStr);
	_delimStr = "";
	getCollapse(); //this will store it's results in the _retStr method.
	_delimStr = oldDelimStr;
	return _retStr;
}

// return a histogram of values and their freqs. in desc. order of frequency
const string &KeyListOpsMethods::getFreqDesc() {
	if (empty()) return _nullVal;

	//for each uniq val, report # occurances, in desc order.
	makeFreqMap();
	//put freq map into multimap where key is the freq and val is the item. In other words, basically a reverse freq map.
	histDescType hist;
	for (; _freqIter != _freqMap.end(); _freqIter++) {
		hist.insert(pair<int, string>(_freqIter->second, _freqIter->first));
	}
	//now iterate through the reverse map we just made and output it's pairs in val:key format.
	_retStr.clear();
	ostringstream s;
	for (histDescType::iterator histIter = hist.begin(); histIter != hist.end(); histIter++) {
		if (histIter != hist.begin()) s << _delimStr;
		s << histIter->second;
		s << ":";
		s << histIter->first;
	}
	_retStr.append(s.str());
	return _retStr;
}
// return a histogram of values and their freqs. in asc. order of frequency
const string &KeyListOpsMethods::getFreqAsc() {
	if (empty()) return _nullVal;

	//for each uniq val, report # occurances, in asc order.
	makeFreqMap();
	//put freq map into multimap where key is the freq and val is the item. In other words, basically a reverse freq map.
	histAscType hist;
	for (; _freqIter != _freqMap.end(); _freqIter++) {
		hist.insert(pair<int, string>(_freqIter->second, _freqIter->first));
//		hist[*(_freqIter->second)] = _freqIter->first;
	}
	//now iterate through the reverse map we just made and output it's pairs in val:key format.
	_retStr.clear();
	ostringstream s;
	for (histAscType::iterator histIter = hist.begin(); histIter != hist.end(); histIter++) {
		if (histIter != hist.begin()) s << _delimStr;
		s << histIter->second;
		s << ":";
		s << histIter->first;
	}
	_retStr.append(s.str());
	return _retStr;
}
// return the first value in the list
const string &KeyListOpsMethods::getFirst() {
	if (empty()) return _nullVal;

	//just the first item.
	begin();
	return getColVal();
}
// return the last value in the list
const string &KeyListOpsMethods::getLast() {
	if (empty()) return _nullVal;

	//just the last item.
	begin();
	for (size_t i = 0; i < getCount() -1; i++) {
		next();
	}
	return getColVal();
}

const string &KeyListOpsMethods::getColVal() {
	const string &retVal = (*_iter)->getField(_column);
	if (_isBam && retVal.empty()) return _nullVal;
	return retVal;
}

double KeyListOpsMethods::getColValNum() {
	const string &strVal = (*_iter)->getField(_column);
	if (!isNumeric(strVal)) {
		_nonNumErrFlag = true;
		ostringstream s;
		_errMsg = " ***** WARNING: Non numeric value ";
		s << strVal;
		s << " in ";
		s << _column;
		s << ".";
		_errMsg.append(s.str());
		return NAN;
	}
	return atof(strVal.c_str());
}

void KeyListOpsMethods::toArray(bool useNum, SORT_TYPE sortVal) {

	//TBD: optimize performance with better memory management.
	if (useNum) {
		_numArray.resize(_keyList->size());
		int i=0;
		for (begin(); !end(); next()) {
			_numArray[i] = getColValNum();
			i++;
		}
	} else {
		_qsArray.resize(_keyList->size());
		int i=0;
		for (begin(); !end(); next()) {
			_qsArray[i] = getColVal();
			i++;
		}
	}
	if (sortVal != UNSORTED) {
		sortArray(useNum, sortVal == ASC);
	}
}

void KeyListOpsMethods::sortArray(bool useNum, bool ascOrder)
{
	if (useNum) {
		if (ascOrder) {
			sort(_numArray.begin(), _numArray.end(), less<double>());
		} else {
			sort(_numArray.begin(), _numArray.end(), greater<double>());
		}
	} else {
		if (ascOrder) {
			sort(_qsArray.begin(), _qsArray.end(), less<string>());
		} else {
			sort(_qsArray.begin(), _qsArray.end(), greater<string>());
		}
	}
}

void KeyListOpsMethods::makeFreqMap() {
	_freqMap.clear();

	//make a map of values to their number of times occuring.
	for (begin(); !end(); next()) {
		_freqMap[getColVal()]++;
	}
	_freqIter = _freqMap.begin();
}
