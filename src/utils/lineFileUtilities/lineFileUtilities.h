#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

using namespace std;

// split a line from a file into a vector of strings.  token = "\t"
//void Tokenize(const string &str, vector<string>& tokens, const string &delimiter = "\t");
//void Tokenize(const string &str, vector<int>& tokens,    const string &delimiter = "\t");

// templated function to convert objects to strings
template <typename T>
inline
std::string ToString(const T & value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

inline
void Tokenize(const string &str, vector<string> &tokens, const string &delimiter = "\t") {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}

inline
void Tokenize(const string &str, vector<int> &tokens, const string &delimiter = "\t") {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiter, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(atoi(str.substr(lastPos, pos - lastPos).c_str()));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}

#endif /* LINEFILEUTILITIES_H */

