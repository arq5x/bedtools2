#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iostream>

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
void Tokenize(const string &str, vector<string> &elems, char delimiter = '\t') 
{
    // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
    // NOTE: this approach intentionally allows consecutive delimiters
    std::stringstream ss(str);
    std::string item;
    while(getline(ss, item, delimiter)) {
        elems.push_back(item);  
    }
}

// tokenize into a list of integers
inline
void Tokenize(const string &str, vector<int> &elems, char delimiter = '\t') 
{

    // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
    // NOTE: this approach intentionally allows consecutive delimiters
    std::stringstream ss(str);
    std::string item;
    while(getline(ss, item, delimiter)) {
        elems.push_back(atoi(item.c_str()));  
    }
}

#endif /* LINEFILEUTILITIES_H */

