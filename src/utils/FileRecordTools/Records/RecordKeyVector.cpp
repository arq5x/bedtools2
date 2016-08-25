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
	_recVec = new vecType();
	_mustDeleteVec = true;
}

RecordKeyVector::RecordKeyVector(Record * item)
: _key(item),
  _currPos(0)
{
	_recVec = new vecType();
	_mustDeleteVec = true;
}

RecordKeyVector::RecordKeyVector(Record * item, vecType *vec)
: _key(item),
  _currPos(0)
{
	_recVec = new vecType(*vec);
	_mustDeleteVec = true;
}

RecordKeyVector::~RecordKeyVector() 
{
	if (_mustDeleteVec) 
	{
		delete _recVec;
		_recVec = NULL;
		_mustDeleteVec = false;
	}
}

const RecordKeyVector &RecordKeyVector::operator=(const RecordKeyVector &other)
{
	setKey(other._key);
	_recVec = other._recVec;
	return *this;
}

RecordKeyVector::iterator_type RecordKeyVector::begin()  {
	_currPos = 0;
	return _recVec->begin();
}

RecordKeyVector::iterator_type RecordKeyVector::next()  {
	_currPos++;
	return _recVec->begin() + _currPos;
}


RecordKeyVector::iterator_type RecordKeyVector::end() {
	return _recVec->end();
}

void RecordKeyVector::erase()  
{
	delete _recVec->at(_currPos);
	_recVec->erase (_recVec->begin() + _currPos);
}

size_t RecordKeyVector::size() const {
	return _recVec->size();
}

bool RecordKeyVector::empty() const {
	return _recVec->empty();
}

void RecordKeyVector::push_back(elemType item) {
	_recVec->push_back(item);
}

Record * RecordKeyVector::getKey() const {
	return _key;
}

void RecordKeyVector::setKey(elemType key) {
	_key = key;
}

void RecordKeyVector::setVector(vecType *vec) {
	_currPos = 0;
	_recVec = vec;
}

void RecordKeyVector::clearVector() {
	_currPos = 0;
	//for (vecType::iterator it = _recVec->begin() ; it != _recVec->end(); ++it)
    //{
    //	delete *it;
    //}
	_recVec->clear();
}

void RecordKeyVector::sortVector() {
	std::sort(_recVec->begin(), _recVec->end(), RecordPtrSortAscFunctor());
}

void RecordKeyVector::swap(RecordKeyVector &other)
{
	//set tmp to this
    elemType tmpKey = _key;
    vecType *tmpVec = _recVec;
    int tmpPos = _currPos;
    bool tmpDel = _mustDeleteVec;

    //set this to other
    _key = other._key;
    _recVec = other._recVec;
    _currPos = other._currPos;
    _mustDeleteVec = other._mustDeleteVec;

    //set other to tmp
    other._key = tmpKey;
    other._recVec = tmpVec;
    other._currPos = tmpPos;
    other._mustDeleteVec = tmpDel;
}


