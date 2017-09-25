/*
 * StrandQueue.h
 *
 *  Created on: Jan 29, 2013
 *      Author: nek3d
 */
#ifndef STRANDQUEUE_H_
#define STRANDQUEUE_H_

#include <vector>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include "Record.h"

using namespace std;

class StrandQueue {
public:
	StrandQueue();
	~StrandQueue();

	Record * top();
	void pop();
	Record * top(Record::strandType strand);
	void pop(Record::strandType strand);
	void push(Record *record);
	size_t size();
	bool empty();

private:
//	static RecordPtrSortFunctor _recSortFunctor;
	typedef priority_queue<Record *, vector<Record *>, RecordPtrSortDescFunctor > queueType;
	vector<queueType *> _queues;
	static const int NUM_QUEUES = 3;

	//we want to be able to iterate over the enumerated strand types in Record.h,
	//which are FORWARD, REVERSE, and UNKNOWN. However, iterating over an enum is hard to
	//do, so we'll use a suggestion found in a forum, and put the enum values into a vector.
	vector<Record::strandType> _strandIdxs;

	int getMinIdx(); //will return the idx of queue with the current min val.

};


#endif // STRANDQUEUE_H_
