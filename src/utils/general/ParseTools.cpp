#include "ParseTools.h"
#include <climits>
#include <cctype>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sstream>


//This functions recognizes only numbers with digits, plus sign, minus sign, decimal point, e, or E. Hexadecimal and pointers not currently supported.
bool isNumeric(const string &str) {
	bool hasDigits = false;
	for (int i=0; i < (int)str.size(); i++) {
		char currChar = str[i];
		if (!(isdigit(currChar) || currChar == '-' || currChar == '.' || currChar == '+' || currChar == 'e' || currChar == 'E')) {
			return false;
		}
		hasDigits |= isdigit(currChar);
	}
	return hasDigits;
}

//As above, but does not allow decimal points
bool isInteger(const string &str) {
	bool hasDigits = false;
	for (int i=0; i < (int)str.size(); i++) {
		char currChar = str[i];
		if (!(isdigit(currChar) || currChar == '-' || currChar == '+' || currChar == 'e' || currChar == 'E')) {
			return false;
		}
		hasDigits |= isdigit(currChar);
	}
	return hasDigits;
}



int str2chrPos(const string &str) {
	return str2chrPos(str.c_str(), str.size());
}

int str2chrPos(const char *str, size_t ulen) {

	if (ulen == 0) {
		ulen = strlen(str);
	}

	//first test for exponents / scientific notation
	for (size_t i=0; i < ulen; i++) {
		if (str[i] == 'e' || str[i] == 'E' || str[i] == '.') {
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

	int sum=0;
	int startPos =0;
	bool isNegative = false;

	//check for negative numbers
	if (str[0] == '-') {
		isNegative = true;
		startPos = 1;
	}

	//check for explicitly positive numbers
	if (str[0] == '+') {
		isNegative = false;
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

bool isHeaderLine(const string &line) {
	if (line[0] == '>') {
		return true;
	}
	if (line[0] == '!') {
		return true;
	}
	if (line[0] == '#') {
		return true;
	}
	
	string tmp = line;
	transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
	//allow chr chrom to start a header line
	if (memcmp(tmp.c_str(), "chrom", 5) == 0 && isspace(tmp[5]) && ! isdigit(tmp[6])) {
		return true;
	}
	//allow chr chrom to start a header line
	if (memcmp(tmp.c_str(), "chr", 3) == 0 && isspace(tmp[3]) && ! isdigit(tmp[4])) {
		return true;
	}
	//UCSC file headers can also start with the words "browser" or "track", followed by a whitespace character.
	if (memcmp(tmp.c_str(), "browser", 7) == 0) {
		return true;
	}
	if (memcmp(tmp.c_str(), "track", 5) == 0) {
		return true;
	}
	if (memcmp(tmp.c_str(), "visibility", 10) == 0) {
		return true;
	}
	return false;
}
