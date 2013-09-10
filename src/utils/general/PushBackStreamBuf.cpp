/*
 * PushBackStream.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: nek3d
 */
#include <cstdio>

#include "PushBackStreamBuf.h"
#include <cstring>
#include <iostream>

PushBackStreamBuf::PushBackStreamBuf(streambuf* primary_stream)
: streambuf(),
  _primary_stream(primary_stream)
 {

}

PushBackStreamBuf::~PushBackStreamBuf()
{
}

int PushBackStreamBuf::sbumpc()
 {
	int c = 0;
	if (!_buffer.empty()) {
		c = _buffer.front()->value();
		_buffer.pop_front();
		return c;
	}
	c = _primary_stream->sbumpc();
	return c;
}

void PushBackStreamBuf::pushBack(const BTlist<int> &newBuf)
{
	_buffer = newBuf;
}

int PushBackStreamBuf::underflow()
{
	int c = 0;

	if(!_buffer.empty()) {
		c=_buffer.front()->value();
		return c;
	}
	c= _primary_stream->sgetc();
	return c;
}
