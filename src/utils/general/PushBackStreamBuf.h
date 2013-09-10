/*
 * PushBackStream.h
 *
 *  Created on: Mar 18, 2013
 *      Author: nek3d
 */

#ifndef PUSHBACKSTREAM_H_
#define PUSHBACKSTREAM_H_

using namespace std;

#include <iostream>
#include "BTlist.h"

class PushBackStreamBuf: public std::streambuf {
public:
	PushBackStreamBuf(streambuf* primary_stream);
	~PushBackStreamBuf();
	void pushBack(const BTlist<int> &vec);

	int sbumpc();
protected:
	int uflow() { return sbumpc(); }

	int underflow();

private:
	streambuf* _primary_stream;
	BTlist<int> _buffer;
	static const int SCAN_BUFFER_SIZE = 4096; //4 KB buffer
	static const int MAIN_BUFFER_SIZE = 128 * 1024; //128K buffer
};


#endif /* PUSHBACKSTREAM_H_ */
