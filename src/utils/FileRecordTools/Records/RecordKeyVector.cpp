/*
 * RecordKeyVector.cpp
 *
 *  Created on: Aug 1, 2014
 *      Author: nek3d
 */


#include "RecordKeyVector.h"
#include <algorithm>

RecordKeyVector::RecordKeyVector()
: _key(NULL),
 _currPos(0)
{
}

RecordKeyVector::RecordKeyVector(const Record * item)
: _key(item),
  _currPos(0)
{
}

RecordKeyVector::RecordKeyVector(const Record * item, const vecType &vec)
: _key(item),
  _currPos(0)
{
	_recVec = vec;
}

RecordKeyVector::~RecordKeyVector() {
}

const RecordKeyVector &RecordKeyVector::operator=(const RecordKeyVector &other)
{
	setKey(other._key);
	_recVec = other._recVec;
	return *this;
}

const RecordKeyVector::const_iterator_type RecordKeyVector::begin()  {
	_currPos = 0;
	return _recVec.begin();
}

const RecordKeyVector::const_iterator_type RecordKeyVector::next()  {
	_currPos++;
	return _recVec.begin() + _currPos;
}


const RecordKeyVector::const_iterator_type RecordKeyVector::end() {
	return _recVec.end();
}

size_t RecordKeyVector::size() const {
	return _recVec.size();
}

bool RecordKeyVector::empty() const {
	return _recVec.empty();
}

void RecordKeyVector::push_back(elemType item) {
	_recVec.push_back(item);
}

const Record *RecordKeyVector::getKey() const {
	return _key;
}

void RecordKeyVector::setKey(elemType key) {
	_key = key;
}

void RecordKeyVector::setVector(const vecType &vec) {
	_recVec = vec;
}

void RecordKeyVector::clearVector() {
	_recVec.clear();
}

void RecordKeyVector::sortVector() {
	std::sort(_recVec.begin(), _recVec.end(), RecordPtrSortAscFunctor());
}


