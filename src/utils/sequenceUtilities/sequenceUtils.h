#ifndef SEQUENCEUTILS_H
#define SEQUENCEUTILS_H

#include <string>
#include <algorithm>
#include <cctype>

using namespace std;

// Performs an in-place sequence reversal
void reverseSequence(string &);

// Performs an in-place reverse complement conversion
void reverseComplement(string &);

// Converts every character in a string to lowercase
void toLowerCase(string &);

// Converts every character in a string to uppercase
void toUpperCase(string &);

#endif
