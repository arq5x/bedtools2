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

bool isNumeric(const string &str);
bool isInteger(const string &str);

//This method is a faster version of atoi, but is limited to a maximum of
//9 digit numbers in Base 10 only. The string may begin with a negative.
//Empty strings, too long strings, or strings containing anything other than
//digits (with the excpetion of a minus sign in the first position)
//will result in error. Errors return INT_MIN.
int str2chrPos(const char *str, size_t len = 0);
int str2chrPos(const string &str);


//int2str is faster but less flexible version of the ToString method in
//lineFileUtilities. Unlike ToString, which uses streams, this method
//can only handle integers. The buffer, which is templated, needs to support the
//assignment operater for char *, meaning it needs a T::operator = (const char *) method.
//strings, strings, stringbuffers, and the like are acceptable.

template<class T>
void int2str(int number, T& buffer, bool appendToBuf = false)
{
	if (number == 0) {
		if (appendToBuf) {
			buffer.append("0");
		} else {
			buffer = "0";
		}
		return;
	}
	//check for negative numbers.
	bool isNegative = number < 0;
	unsigned useNum = number;
	if (isNegative) {
		useNum = 0 - useNum; //convert to positive.
	}

	char tmpBuffer[numeric_limits<int>::digits10 + 3];
	char *tmpBuf = &tmpBuffer[sizeof tmpBuffer];
	*--tmpBuf = '\0';

	while (useNum > 0) {
		*--tmpBuf = (useNum % 10) + '0';
		useNum /= 10;
	}

	if (isNegative) {
		*--tmpBuf = '-';
	}

	if (!appendToBuf) {
		buffer = tmpBuf;
	} else {
		buffer.append(tmpBuf);
	}

}

bool isHeaderLine(const string &line);

string vectorIntToStr(const vector<int> &vec);


//This is a faster version of tokenize that doesn't use strings. Return value is final size of elems vector.
//int Tokenize(const string &str, vector<string> &elems, char delimiter = '\t', int numExpectedItems = 0);

#endif /* PARSETOOLS_H_ */
