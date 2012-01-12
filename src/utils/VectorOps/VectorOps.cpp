/*****************************************************************************
  VectorOps.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Purpose: Container to perform various useful operations 
           (e.g., min, max, count, distinct, etc.) on vectors.

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "VectorOps.h"

// functor to convert a single string to a double
double MakeDouble(const string &element) {
    std::istringstream i(element);
    double x;
    if (!(i >> x)) {
        cerr << "Error: Could not properly convert string to numeric (\"" + element + "\")" << endl;
        exit(1);
    }
    return x;
}

// Constructor
VectorOps::VectorOps(const vector<string> &vec)
: _vecs(vec)
, _vecd()
, _size(vec.size())
{
    _vecd.reserve(vec.size());
}

// Destructor
VectorOps::~VectorOps(void) {
}

double VectorOps::GetSum(void) 
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return accumulate(_vecd.begin(), _vecd.end(), 0.0);
}

double VectorOps::GetMean(void) 
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return accumulate(_vecd.begin(), _vecd.end(), 0.0) / _size;
}

double VectorOps::GetMedian(void) 
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    
    double median = 0.0;
    sort(_vecd.begin(), _vecd.end());
    int totalLines = _vecd.size();
    if ((totalLines % 2) > 0) {
        long mid;
        mid = totalLines / 2;
        median = _vecd[mid];
    }
    else {
        long midLow, midHigh;
        midLow = (totalLines / 2) - 1;
        midHigh = (totalLines / 2);
        median = (_vecd[midLow] + _vecd[midHigh]) / 2.0;
    }
    return median;
}

double VectorOps::GetMin(void)
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return *min_element( _vecd.begin(), _vecd.end() );
}

double VectorOps::GetMax(void)
{
    // convert the vec of strings to a vec of doubles
    transform(_vecs.begin(), _vecs.end(), back_inserter(_vecd), MakeDouble);
    return *max_element( _vecd.begin(), _vecd.end() );
}

uint32_t VectorOps::GetCount(void)
{
    return _vecs.size();
}

uint32_t VectorOps::GetCountDistinct(void)
{
    sort( _vecs.begin(), _vecs.end() );
    _vecs.erase( unique( _vecs.begin(), _vecs.end() ), _vecs.end() );
    return _vecs.size();
}

string VectorOps::GetCollapse(void)
{
    ostringstream collapse;
    for( size_t i = 0; i < _vecs.size(); i++ ) {
        if (i>0)
            collapse << ",";
        collapse << _vecs[i];
    }
    return collapse.str();
}

string VectorOps::GetDistinct(void)
{
    ostringstream distinct;
    // remove duplicate entries from the vector
    // http://stackoverflow.com/questions/1041620/most-efficient-way-to-erase-duplicates-and-sort-a-c-vector
    sort( _vecs.begin(), _vecs.end() );
    _vecs.erase( unique( _vecs.begin(), _vecs.end() ), _vecs.end() );
    
    for( size_t i = 0; i < _vecs.size(); i++ ) {
        if (i>0)
            distinct << ",";
        distinct << _vecs[i];
    }
    return distinct.str();
}