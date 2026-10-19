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
	  _blockSize(blockSize),
	  _growSize(blockSize)
	{
		_freeList.reserve(_blockSize);
	}

	~FreeList() {
		for (int i=0; i < (int)_blocks.size(); i++) {
			delete [] _blocks[i];
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
	vector<T *> _blocks; // one array per growBuffer call
	vector<T *> _freeList;
	int _nextPos;
	int _blockSize;
	int _growSize; // objects in the next block
	// One allocation per block instead of one per object. Nothing deletes a
	// pooled object on its own: deleteObj returns it to _freeList and the
	// destructor frees whole blocks. new T[n]() value-initializes each
	// element, as new T() did. Blocks double from _blockSize up to
	// MAX_BLOCK_MULTIPLE times it, so a pool that keeps growing (database
	// records, never returned) makes few large allocations.
	void growBuffer() {
		if (_growSize <= 0) {
			return;
		}
		T *block = new T[_growSize]();
		_blocks.push_back(block);
		for (int i=0; i < _growSize; ++i) {
			_buffer.push_back(block + i);
		}
		if (_growSize / _blockSize < MAX_BLOCK_MULTIPLE) {
			_growSize *= 2;
		}
	}
	static const int MAX_BLOCK_MULTIPLE = 8;

	//this number determines how many times larger than the blockSize the freeList
	//is allowed to grow. If speed performance is poor, we may need to increase this number.
	//If memory performance is poor, decrease it.
	static const int MAX_MULTIPLE_FREE_LIST_SIZE = 10;
};


#endif /* FREELIST_H_ */
