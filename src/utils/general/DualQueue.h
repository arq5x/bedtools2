/*
 * DualQueue.h
 *
 *  Created on: Jan 29, 2013
 *      Author: nek3d
 */
#ifdef false
#ifndef DUALQUEUE_H_
#define DUALQUEUE_H_

using namespace std;

#include <vector>
#include <queue>
#include <cstdio>
#include <cstdlib>

template <class T> class DualQueueAscending {
public:
	bool operator() ( const T &item1, const T &item2) const {

		printf("\n\nIn comparison method:\n item1=\n");
//		item1->print();
		printf("\nitem2=\n");
//		item2->print();
		printf("\n");

		if( *(item1) < *(item2) ) {
			printf("Item1 less than item2. Returning false.\n");
			return false;
		}
		printf("Item1 not less than item2. Returning true.\n");
		return true;
	}
};

template <class T> class DualQueueDescending {
public:
	bool operator() ( const T &item1, const T &item2) const {
		if( *(item2) < *(item1) ) {
			return false;
		}
		return true;
	}
};


template <class T, template<class T> class CompareFunc> class DualQueue {
public:
	DualQueue() {}
	~DualQueue() {}

	const T & top() const {
		if (empty()) {
			fprintf(stderr, "ERROR. Tried to top from empty dualQueue.\n");
			exit(1);
		}
		if (emptyForward()) {
			return topReverse();
		}
		if (emptyReverse()) {
			return topForward();
		}

		return (topFowardHigherPriorityThanTopReverse() ? topForward() : topReverse());
	}
	void pop() {
		if (empty()) {
			fprintf(stderr, "ERROR. Tried to pop from empty dualQueue.\n");
			exit(1);
		}
		if (emptyForward()) {
			popReverse();
			return;
		}
		if (emptyReverse()) {
			popForward();
			return;
		}
		topFowardHigherPriorityThanTopReverse() ? popForward() : popReverse();
	}
	void push(const T &item) { item->getStrand() ? pushForward(item) : pushReverse(item); }
	size_t size() const { return sizeForward() + sizeReverse(); }
	bool empty() const { return _forwardQueue.empty() && _reverseQueue.empty(); }

	const T & topForward() const { return _forwardQueue.top(); }
	void popForward() { _forwardQueue.pop(); }
	void pushForward(const T &item) { _forwardQueue.push(item); }
	size_t sizeForward() const { return _forwardQueue.size(); }
	bool emptyForward() const { return _forwardQueue.empty(); }


	const T & topReverse() const { return _reverseQueue.top(); }
	void popReverse() { _reverseQueue.pop(); }
	void pushReverse(const T &item) { _reverseQueue.push(item); }
	size_t sizeReverse() const { return _reverseQueue.size(); }
	bool emptyReverse() const { return _reverseQueue.empty(); }


private:
	typedef priority_queue<T, vector<T>, CompareFunc<T> > queueType;
	queueType _forwardQueue;
	queueType _reverseQueue;

	bool topFowardHigherPriorityThanTopReverse() const {
		printf("\n\nIn priority method:\n TopForward=\n");
//		topForward()->print();
		printf("\nTopReverse=\n");
//		topReverse()->print();
		printf("\n");
		if (CompareFunc<T>()(topForward(), topReverse())) {
			printf("Forward higher priority than reverse.\n");
			return true;
		} else {
			printf("Reverse higher priority than forward.\n");
			return false;
		}
	}

};


#endif /* DUALQUEUE_H_ */
#endif
