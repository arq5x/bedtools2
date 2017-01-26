#include "string.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "ParseTools.h"
#include "lineFileUtilities.h"

string::string(size_t capacity)
: _buffer(NULL),
  _currCapacity(capacity),
  _currSize(0)
{
	build();
}

string::string(const string &qs)
:	_buffer(NULL),
	_currCapacity(qs._currCapacity),
	_currSize(0)
{
	build();
	set(qs._buffer, qs._currSize);
}

string::string(const char *inBuf)
{
	size_t len = strlen(inBuf);
	_currCapacity = len +1;

	build();
	set(inBuf, len);
}

string::string(const string &inString)
{
	size_t len = (int)inString.size();
	_currCapacity = len +1;

	build();
	set(inString.c_str(), len);
}

string::string(char c)
{
	_currCapacity =2;

	build();

	char buffer[2];
	buffer[0] = c;
	buffer[1] = 0;

	set(buffer, 1);
}

void string::build() {
	_buffer = (char *)malloc(_currCapacity);
	clear();
}

string::~string(){
	free(_buffer);
}

void string::clear() {
	memset(_buffer, 0, _currCapacity);
	_currSize = 0;
}

void string::release() {
	free(_buffer);
	_currCapacity = DEFAULT_CAPACITY;
	build();
}

string &string::operator = (const char *inBuf){
	set(inBuf, strlen(inBuf));
	return *this;
}

string &string::operator = (const string & inBuf){
	set(inBuf.c_str(), (int)inBuf.size());
	return *this;
}

string &string::operator = (const string & inBuf){
	set(inBuf._buffer, (int)inBuf._currSize);
	return *this;
}

string &string::operator = (char val) {
	clear();
	append(val);
	return *this;
}
string &string::operator = (int val) {
	clear();
	append(val);
	return *this;
}

string &string::operator = (uint32_t val) {
	clear();
	append(val);
	return *this;
}

// string &string::operator = (size_t val) {
// 	clear();
// 	append(val);
// 	return *this;
// }

string &string::operator = (float val) {
	clear();
	append(val);
	return *this;
}

string &string::operator = (double val) {
	clear();
	append(val);
	return *this;
}


string &string::operator += (const string & inBuf)
{
	append(inBuf._buffer, (int)inBuf._currSize);
	return *this;
}

string &string::operator +=(const string &inBuf)
{
	append(inBuf.c_str(), (int)inBuf.size());
	return *this;
}

string &string::operator +=(char c) {

	append(c);
	return *this;
}

string &string::operator += (const char *inBuf)
{
	append(inBuf, strlen(inBuf));
	return *this;
}

string &string::operator += (int num) {
	append(num);
	return *this;
}

string &string::operator += (uint32_t num) {
	append(num);
	return *this;
}

// string &string::operator += (size_t num) {
// 	append(num);
// 	return *this;
// }

string &string::operator += (float num) {
	append(num);
	return *this;
}

string &string::operator += (double num) {
	append(num);
	return *this;
}

bool string::operator == (const string &qs) const {
	if ( _currSize != qs._currSize) {
		return false;
	}
	for (int i= _currSize-1; i > -1; i--) {
		if (_buffer[i] != qs._buffer[i]) return false;
	}
	return true;
}

bool string::operator == (const string &str) const {
	if ( _currSize != str.size()) {
		return false;
	}
	for (int i= _currSize-1; i > -1; i--) {
		if (_buffer[i] != str[i]) return false;
	}
	return true;

}

bool string::stricmp(const string &str) const {
	if (str.size() != _currSize) {
		return true;
	}
	for (size_t i=0; i < _currSize; i++) {
		if (tolower(str[i]) != tolower(_buffer[i])) {
			return true;
		}
	}
	return false;
}

bool string::operator == (const char *str) const {
	size_t inLen = strlen(str);
	if (inLen != _currSize) {
		return false;
	}
	for (int i= _currSize-1; i > -1; i--) {
		if (_buffer[i] != str[i]) return false;
	}
	return true;
}


bool string::operator != (const string &qs) const {
	return !(*this == qs);
}

bool string::operator < (const string &qs) const {
	return (memcmp(_buffer, qs._buffer, max(_currSize, qs._currSize)) < 0);
}

bool string::operator > (const string &qs) const {
	return (memcmp(_buffer, qs._buffer, max(_currSize, qs._currSize))> 0);
}

void string::set(const char *inBuf, size_t newLen) {
	reserve(newLen);
	clear();
	memcpy(_buffer, inBuf, newLen);
	_currSize = newLen;
}

void string::reserve(size_t newLen) {
	newLen++; //always leave room for a null termninator.
	if (_currCapacity <= newLen) {
		while (_currCapacity <= newLen) {
			_currCapacity = _currCapacity << 1;
		}
		_buffer = (char *)realloc(_buffer, _currCapacity );
		if (_buffer == NULL) {
			fprintf(stderr, "Error: failed to reallocate string.\n");
			_currSize = 0;
			_currCapacity = 0;
			exit(1);
		}
		//initialize newly reserved memory.
		memset(_buffer + _currSize, 0, _currCapacity - _currSize);
	}
}

void string::append(char c)
{
	reserve(_currSize +1);
	_buffer[_currSize] = c;
	_currSize++;
}

void string::append(const char *inBuf, size_t inBufLen)
{
	reserve(_currSize + inBufLen);
	memcpy(_buffer + _currSize, inBuf, inBufLen);
	_currSize += inBufLen;
}

void string::append(int num) {
	int2str(num, *this, true);
}

void string::append(uint32_t num) {
 	int2str((int)num, *this, true);
}

// void string::append(size_t num) {
// 	int2str((int)num, *this, true);
// }

void string::append(float num) {
	append(ToString(num));
}

void string::append(double num) {
	append(ToString(num));
}



string &string::assign(const char *inBuf, size_t inBufLen)
{
	clear();
	append(inBuf, inBufLen);
	return *this;
}

void string::resize(size_t newSize, char fillChar)
{
	if (newSize > _currSize) { //grow the string, pad with fillChar
		reserve(newSize);
		memset(_buffer + _currSize, fillChar, newSize -_currSize);
	} else if (newSize < _currSize) { //cut off characters from the end
		memset(_buffer + newSize, 0, _currSize - newSize);
	}
	_currSize = newSize;
}


void string::substr (string &newStr, size_t pos, size_t len) const
{
	if (pos >= _currSize) {
		return;
	}
	if (pos + len >= _currSize) {
		len = _currSize - pos;
	}
	newStr.set(_buffer + pos, len);
}

ostream &operator << (ostream &out, const string &str) {
	out << str._buffer;
	return out;
}
