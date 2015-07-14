#include "ParseTools.h"
#include <climits>
#include <cctype>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sstream>

//This functions recognizes only numbers with digits, plus sign, minus sign, decimal point, e, or E. Hexadecimal and pointers not currently supported.
bool isNumeric(const QuickString &str) {
	for (int i=0; i < (int)str.size(); i++) {
		char currChar = str[i];
		if (!(isdigit(currChar) || currChar == '-' || currChar == '.' || currChar == '+' || currChar == 'e' || currChar == 'E')) {
			return false;
		}
	}
	return true;
}

int str2chrPos(const QuickString &str) {
	return str2chrPos(str.c_str(), str.size());
}

int str2chrPos(const char *str, size_t ulen) {

	if (ulen == 0) {
		ulen = strlen(str);
	}

	//first test for exponents / scientific notation
	bool hasExponent = false;
	for (size_t i=0; i < ulen; i++) {
		if (str[i] == 'e' || str[i] == 'E') {
			std::istringstream ss(str);
			double retVal;
			ss >> retVal;
			return (int)retVal;
		}
	}

	int len=(int)ulen;
	if (len < 1 || len > 10) {
		fprintf(stderr, "***** ERROR: too many digits/characters for integer conversion in string %s. Exiting...\n", str);
		exit(1);
	}

	register int sum=0;
	int startPos =0;
	bool isNegative = false;

	//check for negative numbers
	if (str[0] == '-') {
		isNegative = true;
		startPos = 1;
	}

	for (int i=startPos; i < len; i++) {
		char currChar = str[i];
		if (currChar == 'e' || currChar == 'E') {
			//default to atoi for scientific notation
			return atoi(str);
		}
		if (!isdigit(currChar)) {
			fprintf(stderr, "***** ERROR: illegal character '%c' found in integer conversion of string \"%s\". Exiting...\n", currChar, str);
			exit(1);
		}

		int dig = currChar - 48; //ascii code for zero.
		int power = len -i -1;

		switch (power) {
		case 0:
			sum += dig;
			break;
		case 1:
			sum += dig * 10;
			break;
		case 2:
			sum += dig *100;
			break;
		case 3:
			sum += dig *1000;
			break;
		case 4:
			sum += dig *10000;
			break;
		case 5:
			sum += dig *100000;
			break;
		case 6:
			sum += dig *1000000;
			break;
		case 7:
			sum += dig *10000000;
			break;
		case 8:
			sum += dig *100000000;
			break;
		case 9:
			sum += dig *1000000000;
			break;
		default:
			return 0;
			break;
		}
	}
	return isNegative ? 0 - sum : sum;
}

string vectorIntToStr(const vector<int> &vec) {
	string str;
	str.reserve(vec.size());
	for (int i=0; i < (int)vec.size(); i++) {
		str += (char)(vec[i]);
	}
	return str;
}

bool isHeaderLine(const QuickString &line) {
	if (line[0] == '>') {
		return true;
	}
	if (line[0] == '!') {
		return true;
	}
	if (line[0] == '#') {
		return true;
	}
	//GFF file headers can also start with the words "browser" or "track", followed by a whitespace character.
	if (memcmp(line.c_str(), "browser", 7) == 0 && isspace(line[7])) {
		return true;
	}
	if (memcmp(line.c_str(), "track", 5) == 0 && isspace(line[5])) {
		return true;
	}
	if (memcmp(line.c_str(), "visibility", 10) == 0) {
		return true;
	}
	return false;
}
