/*
 * RecordListPtr.h
 *
 *  Created on: Jun 19, 2014
 *      Author: nek3d
 */

#ifndef RECORD_LIST_H_
#define RECORD_LIST_H_


// Writing our own singly linked list class because:
//
// 1) 	SRecord *L doesn't currently support one that I can find without
// 		C++0x support or SGI extensions.
//
// 2) 	A normal singly linked list doesn't support constant time
// 		insertion or deletion at the end. It also doesn't track
//		your current position or let you delete at that position,
//		but this does.
//
// 3) 	A normal SRecord *L list is doubly-linked, which would do all
// 		these things, but more slowly and with more memory.
//
// 4) 	I can.


#include "FreeList.h"
#include "string.h"
#include <cstring> //for memset
#include "Record.h"

class RecordListNode {
public:
	friend class RecordList;

	RecordListNode() : _next(NULL){}
	RecordListNode(Record * val) : _val(val), _next(NULL) {}
	Record * value() const { return _val; }
	const RecordListNode *next() const { return _next; }
	RecordListNode *next() { return _next; }
	bool hasNext() const { return _next != NULL; }
	void clear() {
		_next = NULL;
		_val = NULL; // WARNING! Since _val is a pointer,
		// the user is responsible for any needed deallocation.
	}

private:

	Record * _val;
	RecordListNode *_next;
};

class RecordList {
public:
	RecordList();

	~RecordList() { clear(); }

	typedef const RecordListNode * const_iterator_type;

	bool empty() const { return _begin == NULL; }
	size_t size() const { return _size; }

	const RecordListNode *front() const { return _begin;}
	RecordListNode *back() const { return _currEnd; }
	void pop_front();
	const RecordListNode *begin() {
		_prevCursor = NULL;
		return _begin;
	}

	const RecordListNode *begin() const { return _begin; }
	const RecordListNode *next();
	const RecordListNode *end() const { return NULL; }

	RecordListNode * deleteCurrent();
	void push_back(Record * &val);

	void clear();
	RecordList &operator=(const RecordList &other);
	void assignNoCopy(RecordList &other);

	void sort() { mergeSort(&_begin); }

private:
	RecordListNode *_begin;
	//keep a cursor to the last element for constant time push_back and pop_back functions.
	RecordListNode *_currEnd;
	RecordListNode *_prevCursor;
	size_t _size;

	bool _dontDelete; //this will usally be false, but set to true when
	//calling assignNoCopy, as our list will actually just be a pointer
	//to another list.



	// The rest of this is just helper methods for sorting.
	// We'll use a mergeSort for optimal time / space efficiency.
	// The following methods are copied almost verbatim from
	// http://www.geeksforgeeks.org/merge-sort-for-linked-list/

	void mergeSort(RecordListNode **headRef);
	RecordListNode* sortedMerge(RecordListNode * a, RecordListNode * b);

	/* Split the nodes of the given list into front and back halves,
	     and return the two lists using the reference parameters.
	     If the length is odd, the extra node should go in the front list.
	     Uses the fast/slow pointer strategy.  */
	void frontBackSplit(RecordListNode* source,
			RecordListNode** frontRef, RecordListNode** backRef);
};


#endif /* RECORD_LIST_H_ */
