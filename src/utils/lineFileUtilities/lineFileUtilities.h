#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>

using namespace std;

// templated function to convert objects to strings
template <typename T>
inline
std::string ToString(const T & value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

// tokenize into a list of strings.
inline
void Tokenize(const string &str, vector<string> &elems, const string &delimiter = "\t") 
{
    char* tok;
    char cchars [str.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, str.c_str());
    tok = strtok(cstr, delimiter.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delimiter.c_str());
    }
}

// tokenize into a list of integers
inline
void Tokenize(const string &str, vector<int> &elems, const string &delimiter = "\t") 
{
    char* tok;
    char cchars [str.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, str.c_str());
    tok = strtok(cstr, delimiter.c_str());
    while (tok != NULL) {
        elems.push_back(atoi(tok));
        tok = strtok(NULL, delimiter.c_str());
    }
}

#endif /* LINEFILEUTILITIES_H */

