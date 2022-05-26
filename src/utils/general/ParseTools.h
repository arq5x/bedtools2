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

using namespace std;

typedef int64_t CHRPOS;

bool isNumeric(const string &str);
bool isInteger(const string &str);

void trimNewlines(string& str);

//This method is a faster version of atoi, but is limited to a maximum of
//9 digit numbers in Base 10 only. The string may begin with a negative.
//Empty strings, too long strings, or strings containing anything other than
//digits (with the excpetion of a minus sign in the first position)
//will result in error. Errors return INT_MIN.
CHRPOS str2chrPos(const char *str, size_t len = 0);
CHRPOS str2chrPos(const string &str);


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
	uint32_t n;
	for(n = 0; number; number /= 10)
		tmp[12 - ++n] = number % 10 + '0';
	if(neg) tmp[12 - ++n] = '-';

	if(appendToBuf) buffer.append(tmp + 12 - n, n);
	else buffer.assign(tmp + 12 - n, n);
}

bool isHeaderLine(const string &line);

string vectorIntToStr(const vector<int> &vec);


//This is a faster version of tokenize that doesn't use strings. Return value is final size of elems vector.
//int Tokenize(const string &str, vector<string> &elems, char delimiter = '\t', int numExpectedItems = 0);

#endif /* PARSETOOLS_H_ */
