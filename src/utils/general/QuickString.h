/*
 * QuickString.h
 *
 *  Created on: Dec 3, 2012
 *      Author: nek3d
 */

#ifndef QUICKSTRING_H_
#define QUICKSTRING_H_

using namespace std;
#include <string>
#include <climits>
#include <ostream>

class QuickString {
public:
	QuickString(size_t capacity = DEFAULT_CAPACITY);
	QuickString(const QuickString &); //need explicit copy constructor
	//so address of _buffer member is not copied! Each QuickString object
	//needs to build it's own buffer.
	QuickString(const char *);
	QuickString(const string &);
	QuickString(char c);
	~QuickString();
	size_t size() const { return _currSize; }
	size_t capacity() const { return _currCapacity; }
	void reserve(size_t newCapacity=0);
	bool empty() const { return _currSize == 0; }

	void clear(); //only clears buffer, doesn't delete it.
	QuickString &operator = (const string &);
	QuickString &operator = (const char *);
	QuickString &operator = (const QuickString &);
	QuickString &operator += (const QuickString &);
	QuickString &operator += (const string &);
	QuickString &operator += (const char *);
	QuickString &operator += (char);

	friend ostream &operator << (ostream &out, const QuickString &str);
	bool operator == (const QuickString &) const;
	bool operator == (const string &) const;
	bool operator == (const char *) const;
	bool operator != (const QuickString &) const;
	bool operator < (const QuickString &) const;
	bool operator > (const QuickString &) const;
	const char *c_str() const { return _buffer; }
	const char &operator [] (int pos) const { return _buffer[pos]; }
	char &operator [] (int pos) { return _buffer[pos]; }

	void append(const QuickString &str) { append(str.c_str(), str.size()); }
	void append(const char *buf, size_t bufLen);
	void append(char c);
	void append(int num);

	QuickString &assign(const char *str, size_t n);
	void resize(size_t n, char c = '\0');
	void substr(QuickString &newStr, size_t pos = 0, size_t len = UINT_MAX) const;

private:
	char *_buffer;
	size_t _currCapacity;
	size_t _currSize;

	static const int DEFAULT_CAPACITY = 256;
	void build();
	void set(const char *len, size_t size);
};


#endif /* QUICKSTRING_H_ */
