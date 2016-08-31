/*
 * string.h
 *
 *  Created on: Dec 3, 2012
 *      Author: nek3d
 */

#ifndef string_H_
#define string_H_

#include <string>
#include <stdint.h>
#include <climits>
#include <ostream>

using namespace std;

class string {
public:
	string(size_t capacity = DEFAULT_CAPACITY);
	string(const string &); //need explicit copy constructor
	//so address of _buffer member is not copied! Each string object
	//needs to build it's own buffer.
	string(const char *);
	string(const string &);
	string(char c);
	~string();
	size_t size() const { return _currSize; }
	size_t capacity() const { return _currCapacity; }
	void reserve(size_t newCapacity=0);
	bool empty() const { return _currSize == 0; }

	void clear(); //only clears buffer, doesn't delete it.
	void release(); //will deallocate current buffer, reallocate it at default size.
	string &operator = (const string &);
	string &operator = (const char *);
	string &operator = (const string &);
	string &operator = (char);
	string &operator = (int);
	string &operator = (uint32_t);
	//string &operator = (size_t);
	string &operator = (float);
	string &operator = (double);
	string &operator += (const string &);
	string &operator += (const string &);
	string &operator += (const char *);
	string &operator += (char);
	string &operator += (int);
	string &operator += (uint32_t);
	//string &operator += (size_t);
	string &operator += (float);
	string &operator += (double);

	friend ostream &operator << (ostream &out, const string &str);
	bool operator == (const string &) const;
	bool operator == (const string &) const;
	bool operator == (const char *) const;
	bool operator != (const string &) const;
	bool operator < (const string &) const;
	bool operator > (const string &) const;
	const char *c_str() const { return _buffer; }
	const string str() const { return _buffer; }
	const char &operator [] (int pos) const { return _buffer[pos]; }
	char &operator [] (int pos) { return _buffer[pos]; }
	char &at(size_t pos) { return _buffer[pos]; }
	bool stricmp(const string &str) const; //case insensitive compare. False if same aside from case, true if different aside from case.

	void append(const string &str) { append(str.c_str(), str.size()); }
	void append(const char *buf, size_t bufLen);
	void append(char c);

	//These are not templated because float and double require a stringstream based
	//implementation, while the integer append uses a much faster home-brewed algorithm
	//for better performance.
	void append(int num);
	void append(uint32_t num);
	//void append(size_t num);
	void append(float num);
	void append(double num);



	string &assign(const char *str, size_t n);
	void resize(size_t n, char c = '\0');
	void substr(string &newStr, size_t pos = 0, size_t len = UINT_MAX) const;

private:
	char *_buffer;
	size_t _currCapacity;
	size_t _currSize;

	static const int DEFAULT_CAPACITY = 8;
	void build();
	void set(const char *len, size_t size);
};


#endif /* string_H_ */
