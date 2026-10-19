/*
 * parseTools.h
 *
 *  Created on: Dec 3, 2012
 *      Author: nek3d
 */

#ifndef PARSETOOLS_H_
#define PARSETOOLS_H_

#include <cstring> //for memset
#include <limits>
#include <string>
#include <algorithm>
#include <vector>
#include "string.h"
#include <cstdio>
#include <cstdlib>
#include "BedtoolsTypes.h"

using namespace std;

bool isNumeric(const string &str);
bool isInteger(const string &str);

void trimNewlines(string& str);

//This method is a faster version of atoi, but is limited to a maximum of
//9 digit numbers in Base 10 only. The string may begin with a negative.
//Empty strings, too long strings, or strings containing anything other than
//digits (with the excpetion of a minus sign in the first position)
//will result in error. Errors return INT_MIN.
CHRPOS str2chrPos(const char *str, size_t len = 0);

// Whether the 8 bytes at s are ASCII digits; if so their value goes to v.
// The bytes are assembled in text order whatever the byte order (compilers
// make it one load on little-endian targets), then converted in-word: digit
// pairs, then groups of 4, then of 8.
inline bool digits8(const char *s, uint64_t &v) {
	const unsigned char *u = (const unsigned char *)s;
	uint64_t x = (uint64_t)u[0] | (uint64_t)u[1] << 8 | (uint64_t)u[2] << 16 | (uint64_t)u[3] << 24
		| (uint64_t)u[4] << 32 | (uint64_t)u[5] << 40 | (uint64_t)u[6] << 48 | (uint64_t)u[7] << 56;
	const uint64_t high = 0xF0F0F0F0F0F0F0F0ULL, zeros = 0x3030303030303030ULL;
	// every byte in 0x30-0x3F, and still so after adding 6 (no carry once
	// the first test holds), is a byte in 0x30-0x39
	if ((x & high) != zeros || ((x + 0x0606060606060606ULL) & high) != zeros) {
		return false;
	}
	x -= zeros;
	x = (x * 10 + (x >> 8)) & 0x00FF00FF00FF00FFULL;
	x = (x * 100 + (x >> 16)) & 0x0000FFFF0000FFFFULL;
	x = (x * 10000 + (x >> 32)) & 0xFFFFFFFFULL;
	v = x;
	return true;
}

// Whether s[0, n) is 1 to 18 ASCII digits; if so its value, the one
// str2chrPos gives it, goes to val. Anything else (empty, a sign, an
// exponent, a NUL, 19 or more digits) is left to str2chrPos and its error.
inline bool str2chrPosDigits(const char *s, size_t n, CHRPOS &val) {
	if (n - 1 >= 18) {
		return false;
	}
	uint64_t v = 0, w;
	size_t i = 0;
	if (n >= 8) {
		if (!digits8(s, v)) {
			return false;
		}
		i = 8;
		if (n >= 16) {
			if (!digits8(s + 8, w)) {
				return false;
			}
			v = v * 100000000ULL + w;
			i = 16;
		}
	}
	for (; i < n; i++) {
		unsigned d = (unsigned char)s[i] - (unsigned)'0';
		if (d > 9) {
			return false;
		}
		v = v * 10 + d;
	}
	val = (CHRPOS)v;
	return true;
}

inline CHRPOS str2chrPos(const string &str) {
	CHRPOS val;
	if (str2chrPosDigits(str.data(), str.size(), val)) {
		return val;
	}
	return str2chrPos(str.c_str(), str.size());
}


//int2str is faster but less flexible version of the ToString method in
//lineFileUtilities. Unlike ToString, which uses streams, this method
//can only handle integers. The buffer, which is templated, needs to support the
//assignment operater for char *, meaning it needs a T::operator = (const char *) method.
//strings, strings, stringbuffers, and the like are acceptable.

template<class T, class U>
void int2str(U number, T& buffer, bool appendToBuf = false)
{
	if(number == 0)
	{
		if(appendToBuf) buffer.append("0");
		else buffer.assign("0", 1);
		return;
	}
	char tmp[12];

	bool neg = number < 0;
	if(neg) number = -number;
	unsigned int n;
	for(n = 0; number; number /= 10)
		tmp[12 - ++n] = number % 10 + '0';
	if(neg) tmp[12 - ++n] = '-';

	if(appendToBuf) buffer.append(tmp + 12 - n, n);
	else buffer.assign(tmp + 12 - n, n);
}

bool isHeaderLine(const string &line);
// As above for the len bytes at line, which must be followed by a NUL.
bool isHeaderLine(const char *line, size_t len);

string vectorIntToStr(const vector<int> &vec);


//This is a faster version of tokenize that doesn't use strings. Return value is final size of elems vector.
//int Tokenize(const string &str, vector<string> &elems, char delimiter = '\t', int numExpectedItems = 0);

#endif /* PARSETOOLS_H_ */
