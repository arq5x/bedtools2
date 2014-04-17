/*
 * parseTools.h
 *
 *  Created on: Dec 3, 2012
 *      Author: nek3d
 */

#ifndef PARSETOOLS_H_
#define PARSETOOLS_H_

using namespace std;

#include <cstring> //for memset
#include <string>
#include <vector>
#include "QuickString.h"

bool isNumeric(const QuickString &str);

//This method is a faster version of atoi, but is limited to a maximum of
//9 digit numbers in Base 10 only. The string may begin with a negative.
//Empty strings, too long strings, or strings containing anything other than
//digits (with the excpetion of a minus sign in the first position)
//will result in error. Errors return INT_MIN.
int str2chrPos(const char *str, size_t len = 0);
int str2chrPos(const QuickString &str);


//int2str is faster but less flexible version of the ToString method in
//lineFileUtilities. Unlike ToString, which uses streams, this method
//can only handle integers. The buffer, which is templated, needs to support the
//assignment operater for char *, meaning it needs a T::operator = (const char *) method.
//strings, QuickStrings, stringbuffers, and the like are acceptable.

template<class T>
void int2str(int number, T& buffer, bool appendToBuf = false)
{

	register int useNum = number;
	if (useNum == 0) {
		if (appendToBuf) {
			buffer.append("0");
		} else {
			buffer = "0";
		}
		return;
	}
	//check for negative numbers.
	bool isNegative = useNum < 0;
	if (isNegative) {
		useNum = 0 - useNum; //convert to positive.
	}

	//figure out how many digits we have
	register int power = 10;
	int numChars = 2 + (isNegative ? 1: 0);
	while (power  <= useNum) {
		power *= 10;
		numChars++;
	}

	char tmpBuf[numChars];
	memset(tmpBuf, 0, numChars);

	int bufIdx=0;
	if (isNegative) {
		tmpBuf[0] = '-';
		bufIdx = 1;
	}
	register int currDig=0;

	power /= 10;
	while (power) {

		currDig = useNum / power;
		useNum -= currDig * power;
		tmpBuf[bufIdx] = currDig + 48; //48 is ascii for zero.
		bufIdx++;
		power /= 10;
	}
	if (!appendToBuf) {
		buffer = tmpBuf;
	} else {
		buffer.append(tmpBuf, numChars-1);
	}

}

bool isHeaderLine(const QuickString &line);

string vectorIntToStr(const vector<int> &vec);


//This is a faster version of tokenize that doesn't use strings. Return value is final size of elems vector.
int Tokenize(const QuickString &str, vector<QuickString> &elems, char delimiter = '\t', int numExpectedItems = 0);

#endif /* PARSETOOLS_H_ */
