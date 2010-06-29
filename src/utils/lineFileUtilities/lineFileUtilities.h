#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <algorithm>


using namespace std;

// split a line from a file into a vector of strings.  token = "\t"
void Tokenize(const string &str, vector<string>& tokens, const string &delimiter = "\t");
void Tokenize(const string &str, vector<int>& tokens,    const string &delimiter = "\t");


#endif /* LINEFILEUTILITIES_H */
