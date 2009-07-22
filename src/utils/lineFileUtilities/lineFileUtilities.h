#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <algorithm>


using namespace std;

// split a line from a file into a vector of strings.  token = "\t"
void Tokenize(const string& str, vector<string>& tokens);


#endif /* LINEFILEUTILITIES_H */
