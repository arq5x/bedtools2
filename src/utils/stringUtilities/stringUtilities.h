#ifndef STRINGUTILITIES_H
#define STRINGUTILITIES_H

#include <iostream>
#include <string>
#include <cctype>
#include <cwctype>
#include <stdexcept>

using namespace std;

/****************************************************
// isInteger(s): Tests if string s is a valid integer
*****************************************************/
inline bool isInteger(const std::string& s) {
    int len = s.length();
    for (int i = 0; i < len; i++)
        if (!std::isdigit(s[i])) return false;
    return true;
}

string toLower(const string& s) {
	string t = s;
	for (string::iterator p = t.begin(); p != t.end(); ++p) 
	{
		*p = tolower(*p);
	}
	return t;
}

#endif /* STRINGUTILITIES_H */

