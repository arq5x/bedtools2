#ifndef STRINGUTILITIES_H
#define STRINGUTILITIES_H

#include <cctype>
#include <string>

/****************************************************
// isInteger(s): Tests if string s is a valid integer
*****************************************************/
inline bool isInteger(const std::string& s) {
    int len = s.length();
    for (int i = 0; i < len; i++) {
        if (!std::isdigit(s[i])) return false;
    return true;
}

#endif /* STRINGUTILITIES_H */

