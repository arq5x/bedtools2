/*
 * PushBackStream.h
 *
 *  Created on: Mar 18, 2013
 *      Author: nek3d
 */

#ifndef PUSHBACKSTREAM_H_
#define PUSHBACKSTREAM_H_

#include <iostream>
#include "BTlist.h"

using namespace std;

class PushBackStreamBuf: public std::streambuf {
public:
	friend class InputStreamMgr;
	PushBackStreamBuf(streambuf* primary_stream);
	~PushBackStreamBuf();
	void pushBack(const BTlist<int> &vec);
//	void push_front(const BTlist<int> &vec) { _buffer.push_front(vec); }

	int sbumpc();
	void clear() { _buffer.clear(); }
protected:
	int uflow() { return sbumpc(); }

	int underflow();

private:
	streambuf* _primary_stream;
	BTlist<int> _buffer;
};


#endif /* PUSHBACKSTREAM_H_ */
