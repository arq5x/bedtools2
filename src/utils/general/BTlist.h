/*
 * BTlist.h
 *
 *  Created on: Jan 18, 2013
 *      Author: nek3d
 */

#ifndef BT_LIST_H_
#define BT_LIST_H_

//
// Writing our own singly linked list template class because:
//
// 1) 	STL doesn't currently support one that I can find without
// 		C++0x support or SGI extensions.
//
// 2) 	A normal singly linked list doesn't support constant time
// 		insertion or deletion at the end. It also doesn't track
//		your current position or let you delete at that position,
//		but this does.
//
// 3) 	A normal STL list is doubly-linked, which would do all
// 		these things, but more slowly and with more memory.
//
// 4) 	I can.


#include "FreeList.h"
#include "string.h"
#include <sstream>
#include <stack>
#include <cstring> //for memset

template <class T> class BTlist;

template <class T>
class BTlistNode {
public:
	friend class BTlist<T>;

	BTlistNode() : _next(NULL){}
	BTlistNode(const T &val) : _val(val), _next(NULL) {}
	const T &value() const { return _val; }
	const BTlistNode<T> *next() const { return _next; }
	BTlistNode<T> *next() { return _next; }
	bool hasNext() const { return _next != NULL; }
	void clear() {
		_next = NULL;
		//for the value, set numbers to 0, ptrs to NULL, etc.
		memset((void *)_val, 0, sizeof(T));
	}

private:

	T _val;
	BTlistNode *_next;
};

template <class T> class BTlist {
public:
	BTlist() :
		_begin(NULL),
		_currEnd(NULL),
		_prevCursor(NULL),
		_size(0),
		_dontDelete(false)
	{
	}
	//This version of the constructor will convert a string into a BTlist.
	//It is effectively the opposite of the toStr method.
	BTlist(const string &str) :
		_begin(NULL),
		_currEnd(NULL),
		_prevCursor(NULL),
		_size(0),
		_dontDelete(false)
	{
		for (int i=0; i < (int)str.size(); i++) {
			push_back((T)(str[i]));
		}
	}

	~BTlist() {
		clear();
	}

	typedef const BTlistNode<T> * const_iterator_type;

	bool empty() const { return _begin == NULL; }
	size_t size() const { return _size; }

	const BTlistNode<T> *front() const { return _begin;}
	BTlistNode<T> *back() const { return _currEnd; }
	void pop_front() {
		if (empty()) {
			return;
		}
		BTlistNode<T> *killNode = _begin;
		if (_begin->_next == NULL) { //this is the only item. List is 1.
			delete killNode;
			_begin = NULL;
			_currEnd = NULL;
			_prevCursor = NULL;
			_size = 0;
			return;
		}
		if (_prevCursor == _begin) {
			_prevCursor = _begin->_next;
		}
		_begin = _begin->next();
		delete killNode;
		_size--;
	}

	const BTlistNode<T> *begin() {
		_prevCursor = NULL;
		return _begin;
	}

	const BTlistNode<T> *begin() const {
		return _begin;
	}
	const BTlistNode<T> *next() {
		if (_prevCursor == NULL) {
			_prevCursor = _begin;
			return _begin->_next;
		}
		_prevCursor = _prevCursor->_next;
		return _prevCursor->_next;
	}
	const BTlistNode<T> *end() const { return NULL; }

	BTlistNode<T> * deleteCurrent() {
		//delete what the current cursor is pointing to,
		//meaning whatever the prevCursor's next member points to.
		BTlistNode<T> *returnNode = NULL;
		BTlistNode<T> *killNode = NULL;
		if (_prevCursor == NULL) {
			//deleting first item in list.
			killNode = _begin;
			_begin = _begin->_next;
			returnNode = _begin;

		} else {
			killNode = _prevCursor->_next;
			_prevCursor->_next = _prevCursor->_next->_next;
			returnNode = _prevCursor->_next;
		}
		if (_currEnd == killNode) {
			//deleting last item in list
			_currEnd = _prevCursor; //back up the current end.
		}
		delete killNode;
		_size--;

		return returnNode;
	}

	void push_back(const T &val) {
		BTlistNode<T> *newNode = new BTlistNode<T>(val);
		if (empty()) {
			_begin = newNode;
			_currEnd = newNode;
			_prevCursor = NULL;
		} else {
			_currEnd->_next = newNode;
			_prevCursor = _currEnd;
			_currEnd = newNode;
		}
		_size++;
	}

//	void push_front(const T &val) {
//		BTlistNode<T> *newNode = new BTlistNode<T>(val);
//		if (empty()) {
//			_begin = newNode;
//			_currEnd = newNode;
//			_prevCursor = NULL;
//		} else {
//			newNode->_next = _begin;
//			_begin = newNode;
//			_prevCursor = _currEnd;
//		}
//		_size++;
//	}
//	void push_front(BTlist<T> &newList) {
//		stack<T> theStack;
//		for (const_iterator_type iter = newList.begin(); iter != newList.end(); iter = newList.next()) {
//			theStack.push(iter->_val);
//		}
//		while (!theStack.empty()) {
//			const T &val = theStack.top();
//			theStack.pop();
//			push_front(val);
//		}
//	}

	void clear() {
		if (_dontDelete) {
			return;
		}
		//delete all nodes, set cursors to NULL, set size to 0.
		BTlistNode<T> *killNode = _begin;
		while (killNode != NULL) {
			BTlistNode<T> *nextNode = killNode->_next;
			delete killNode;
			killNode = nextNode;
		}
		_begin = NULL;
		_currEnd = NULL;
		_prevCursor = NULL;
		_size = 0;
	}

	const BTlist<T> &operator=(const BTlist<T> &other) {
		//need to make copies of all nodes and assign values
		clear();

		BTlistNode<T> *node = other._begin;
		while (node != NULL) {
			push_back(node->_val);
			node = node->_next;
		}
		return *this;
	}
	void assignNoCopy(BTlist<T> &other) {
		_begin = other._begin;
		_size = other._size;
		_currEnd = other._currEnd;
		_prevCursor = other._prevCursor;
		_dontDelete = true;
	}

	//this method will convert the contents of the list into a string.
	//It assumes that the templated type of the listNode can be converted to a
	//char using the (char) cast. The user must ensure that this is the case.
	//The 2nd parameter needs only be true if you wish to append the data
	//to the string's existing contents. Otherwise, it will be cleared first.
	void toStr(string &str, bool append = false) const 
	{
		if (!append) {
			str.clear();
		}
		str.reserve(_size + str.size());
		for (const BTlistNode<T> *iter = begin(); iter != end(); iter = iter->next()) {
			str += (char)iter->value();
		}
	}

private:
	BTlistNode<T> *_begin;
	//keep a cursor to the last element for constant time push_back and pop_back functions.
	BTlistNode<T> *_currEnd;
	BTlistNode<T> *_prevCursor;
	size_t _size;

	bool _dontDelete; //this will usally be false, but set to true when
	//calling assignNoCopy, as our list will actually just be a pointer
	//to another list.


};

#endif /* BT_LIST_H_ */
