/*
 * FreeList.h
 *
 *  Created on: Nov 30, 2012
 *      Author: nek3d
 */

#ifndef FREELIST_H_
#define FREELIST_H_

#include <cstddef> //defines NULL
#include <deque>
#include <vector>

using namespace std;

template <class T>
class FreeList {
public:
	FreeList(int blockSize=512)
	: _nextPos(0),
	  _blockSize(blockSize)
	{
		_freeList.reserve(_blockSize);
	}

	~FreeList() {
		for (int i=0; i < (int)_buffer.size(); i++) {
			delete _buffer[i];
		}
	}

	void clear() {
		_nextPos = 0;
		_freeList.clear();
	}

	T *newObj() {
		T *ptr = NULL;

		if (_freeList.size() > 0) {
			ptr = _freeList.back();
			_freeList.pop_back();
		} else {
			if (_nextPos == (int)_buffer.size()) {
				growBuffer();
			}
			ptr = _buffer[_nextPos];
			_nextPos ++;
		}

		return ptr;
	}

	void deleteObj(const T *ptr) {
		if (ptr != NULL) {
			T *volatilePtr = const_cast<T *>(ptr); //remove const-ness.
			volatilePtr->clear();
			_freeList.push_back(volatilePtr);
		}
	}

	int capacity() const { return _buffer.size(); }

private:
	deque<T *> _buffer;
	vector<T *> _freeList;
	int _nextPos;
	int _blockSize;
	void growBuffer() {
		for (int i=0; i < _blockSize; ++i) {
			_buffer.push_back(new T());
		}
	}

	//this number determines how many times larger than the blockSize the freeList
	//is allowed to grow. If speed performance is poor, we may need to increase this number.
	//If memory performance is poor, decrease it.
	static const int MAX_MULTIPLE_FREE_LIST_SIZE = 10;
};


#endif /* FREELIST_H_ */
