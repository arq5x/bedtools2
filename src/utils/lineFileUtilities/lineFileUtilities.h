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

// tokenize a column string into a list of integers.
inline
void TokenizeColumns(const string &str, vector<int> &elems) 
{
    vector<string> col_sets;
    Tokenize(str, col_sets, ',');

    for( size_t i = 0; i < col_sets.size(); i++ ) {
        string col_set = col_sets[i];
        if( string::npos == col_set.find("-") ){
            elems.push_back(atoi(col_set.c_str()));
        }
        else {
            vector<string> ends;
            Tokenize(col_set, ends, '-');
            int start = atoi(ends[0].c_str());
            int end = atoi(ends[1].c_str());
            if(start <= end){
                for(int i = start; i <= end; i++){
                    elems.push_back(i);
                }
            }
            else {
                for(int i = start; i >= end; i--){
                    elems.push_back(i);
                }
            } 
        }
    }
}


#endif /* LINEFILEUTILITIES_H */

