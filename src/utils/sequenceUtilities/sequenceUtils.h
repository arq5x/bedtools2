#ifndef SEQUENCEUTILS_H
#define SEQUENCEUTILS_H

#include <string>
#include <algorithm>
#include <cctype>

using namespace std;

// Performs an in-place sequence reversal
void reverseSequence(string &seq);

// Performs an in-place reverse complement conversion
void reverseComplement(string &seq);

// Converts every character in a string to lowercase
void toLowerCase(string &seq);

// Converts every character in a string to uppercase
void toUpperCase(string &seq);

// Calculates the number of a, c, g, t, n, and other bases found in a sequence
void getDnaContent(const string &seq, int &a, int &c, int &g, int &t, int &n, int &other);

int countPattern(const string &seq, const string &pattern, bool ignoreCase);

#endif
