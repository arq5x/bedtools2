/*
 * StrandQueue.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nek3d
 */

#include "StrandQueue.h"

StrandQueue::StrandQueue() {

	for (int i=0; i < NUM_QUEUES; i++) {
		queueType *queue = new queueType();
		_queues.push_back(queue);
	}
	_strandIdxs.resize(3);
	_strandIdxs[0] = Record::FORWARD;
	_strandIdxs[1] = Record::REVERSE;
	_strandIdxs[2] = Record::UNKNOWN;
}

StrandQueue::~StrandQueue() {
	for (int i=0; i < NUM_QUEUES; i++) {
		delete _queues[i];
	}
}

Record *StrandQueue::top()
{
	int minIdx = getMinIdx();
	if (minIdx == -1) return NULL;
	return _queues[minIdx]->top();
}

void StrandQueue::pop() {
	int minIdx = getMinIdx();
	if (minIdx == -1) return;
	_queues[minIdx]->pop();
}

Record * StrandQueue::top(Record::strandType strand) {
	Record *record = NULL;
	switch (strand) {
	case Record::FORWARD:
		if (_queues[0]->empty()) return NULL;
		record = _queues[0]->top();
		break;
	case Record::REVERSE:
		if (_queues[1]->empty()) return NULL;

		record = _queues[1]->top();
		break;
	case Record::UNKNOWN:
		if (_queues[0]->empty()) return NULL;
		record = _queues[2]->top();
		break;
	default:
		break;
	}
	return record;
}

void StrandQueue::pop(Record::strandType strand) {
	switch (strand) {
	case Record::FORWARD:
		if (_queues[0]->empty()) return;
		_queues[0]->pop();
		break;
	case Record::REVERSE:
		if (_queues[1]->empty()) return;
		_queues[1]->pop();
		break;
	case Record::UNKNOWN:
		if (_queues[2]->empty()) return;
		_queues[2]->pop();
		break;
	default:
		break;
	}
}

void StrandQueue::push(Record *record) {
	switch (record->getStrandVal()) {
	case Record::FORWARD:
		_queues[0]->push(record);
		break;
	case Record::REVERSE:
		_queues[1]->push(record);
		break;
	case Record::UNKNOWN:
		_queues[2]->push(record);
		break;
	default:
		break;
	}
}

size_t StrandQueue::size() {
	size_t sumSize = 0;
	for (int i = 0; i < NUM_QUEUES; i++) {
		sumSize += _queues[i]->size();
	}
	return sumSize;
}

bool StrandQueue::empty() {
	for (int i = 0; i < NUM_QUEUES; i++) {
		if (!_queues[i]->empty()) {
			return false;
		}
	}
	return true;
}


int StrandQueue::getMinIdx() {
	if (empty()) return -1;
	const Record *minRec = NULL;
	int minIdx = -1;
	for (int i = 0; i < NUM_QUEUES; i++) {
		if (_queues[i]->empty()) continue;
		const Record *currTop = _queues[i]->top();
		if (currTop == NULL) continue;
		if (minRec == NULL || *currTop < *minRec) {
			minRec = currTop;
			minIdx = i;
		}
	}
	return minIdx;
}

