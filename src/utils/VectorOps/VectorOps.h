/*****************************************************************************
  VecOps.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef VECTOROPS_H
#define VECTOROPS_H

#include <vector>
#include <map>
#include <math.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits.h>

using namespace std;


//************************************************
// VectorOps Class methods and elements
//************************************************
class VectorOps {

public:

    // Constructor. Initialize with a vector of strings.
    VectorOps(const vector<string> &vec);

    // Destructor
    ~VectorOps(void);

    // user-interface
    double GetSum(void);              // return the total of the values in the vector
    double GetMean(void);             // return the average value in the vector
    double GetStddev(void);           // return the standard deviation
    double GetSstddev(void);          // return the sample standard deviation 
    double GetMedian(void);           // return the median value in the vector
    string GetMode(void);             // return the most common value in the vector
    string GetAntiMode(void);         // return the least common value in the vector
    double GetMin(void);              // return the minimum element of the vector
    double GetMax(void);              // return the maximum element of the vector
    uint32_t GetCount(void);          // return the count of element in the vector
    uint32_t GetCountDistinct(void);  // return a the count of _unique_ elements in the vector
    string GetCollapse(void);         // return a comma-separated list of elements
    string GetConcat(void);           // return a concatenation of all elements in the vector
    string GetDistinct(void);         // return a comma-separated list of the _unique_ elements
    string GetFreqDesc(void);         // return a histogram of values and their frequencies in desc. order of frequency
    string GetFreqAsc(void);          // return a histogram of values and their frequencies in asc. order of frequency

private:
    vector<string> _vecs;
    vector<double> _vecd;
    uint32_t _size;
};

#endif /* VECTOROPS_H */
