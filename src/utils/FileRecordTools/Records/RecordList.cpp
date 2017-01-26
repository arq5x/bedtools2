/*
 * RecordList.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: nek3d
 */

#include "RecordList.h"

RecordList::RecordList() :
		_begin(NULL),
		_currEnd(NULL),
		_prevCursor(NULL),
		_size(0),
		_dontDelete(false)
{
}

void RecordList::pop_front() {
	if (empty()) {
		return;
	}
	RecordListNode *killNode = _begin;
	if (_begin->_next == NULL) { //this is the only item. List size is 1.
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

const RecordListNode *RecordList::next() {
	if (_prevCursor == NULL) {
		_prevCursor = _begin;
		return _begin->_next;
	}
	_prevCursor = _prevCursor->_next;
	return _prevCursor->_next;
}


RecordListNode *RecordList::deleteCurrent() {
	//delete what the current cursor is pointing to,
	//meaning whatever the prevCursor's next member points to.
	RecordListNode *returnNode = NULL;
	RecordListNode *killNode = NULL;
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

void RecordList::push_back(Record * &val) {
	RecordListNode *newNode = new RecordListNode(val);
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


void RecordList::clear() {
        _size = 0;
	if (_dontDelete) {
		return;
	}
	//delete all nodes, set cursors to NULL, set size to 0.
	RecordListNode *killNode = _begin;
	while (killNode != NULL) {
		RecordListNode *nextNode = killNode->_next;
		delete killNode;
		killNode = nextNode;
	}
	_begin = NULL;
	_currEnd = NULL;
	_prevCursor = NULL;
	_size = 0;
}

RecordList &RecordList::operator=(const RecordList &other) {
	//need to make copies of all nodes and assign values
	clear();

	RecordListNode *node = other._begin;
	while (node != NULL) {
		push_back(node->_val);
		node = node->_next;
	}
	return *this;
}
void RecordList::assignNoCopy(RecordList &other) {
	_begin = other._begin;
	_size = other._size;
	_currEnd = other._currEnd;
	_prevCursor = other._prevCursor;
	_dontDelete = true;
}


void RecordList::mergeSort(RecordListNode **headRef) {
	RecordListNode *head = *headRef;
	RecordListNode *a;
	RecordListNode *b;

	/* Base case -- length 0 or 1 */
	if ((head == NULL) || (head->_next == NULL)) {
		return;
	}

	/* Split head into 'a' and 'b' sublists */
	frontBackSplit(head, &a, &b);

	/* Recursively sort the sublists */
	mergeSort(&a);
	mergeSort(&b);

	/* answer = merge the two sorted lists together */
	*headRef = sortedMerge(a, b);
}

RecordListNode* RecordList::sortedMerge(RecordListNode * a, RecordListNode * b) {
	RecordListNode *result = NULL;

	/* Base cases */
	if (a == NULL) return(b);
	else if (b==NULL) return(a);

	/* Pick either a or b, and recur */
	if (a->_val->lessThan(b->_val)) {
		result = a;
		result->_next = sortedMerge(a->_next, b);
	} else {
		result = b;
		result->_next = sortedMerge(a, b->_next);
	}
	return(result);
}

/* URecord *ILIRecord *Y FUNCRecord *IONS */
/* Split the nodes of the given list into front and back halves,
	 and return the two lists using the reference parameters.
	 If the length is odd, the extra node should go in the front list.
	 Uses the fast/slow pointer strategy.  */
void RecordList::frontBackSplit(RecordListNode* source,
		RecordListNode** frontRef, RecordListNode** backRef) {
	RecordListNode *fast;
	RecordListNode *slow;
	if (source==NULL || source->_next==NULL) {
		/* length < 2 cases */
		*frontRef = source;
		*backRef = NULL;
	} else {
		slow = source;
		fast = source->_next;

		/* Advance 'fast' two nodes, and advance 'slow' one node */
		while (fast != NULL) {
			fast = fast->_next;
			if (fast != NULL) {
				slow = slow->_next;
				fast = fast->_next;
			}
		}

		/* 'slow' is before the midpoint in the list, so split it in two
		at that point. */
		*frontRef = source;
		*backRef = slow->_next;
		slow->_next = NULL;
	}
}



