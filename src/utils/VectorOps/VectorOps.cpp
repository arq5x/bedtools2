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

struct ValueGreaterThan
{
    bool operator()( const vector< pair<int, string> >::value_type& lhs,
        const vector< pair<int, string> >::value_type& rhs ) const
    {
        return lhs.first > rhs.first;
    }
};

struct ValueLessThan
{
    bool operator()( const vector< pair<int, string> >::value_type& lhs,
        const vector< pair<int, string> >::value_type& rhs ) const
    {
        return lhs.first < rhs.first;
    }
};

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

double VectorOps::GetStddev(void) 
{
    double mean = GetMean();
    // get the variance
    double totalVariance = 0.0;
    vector<double>::const_iterator dIt  = _vecd.begin();
    vector<double>::const_iterator dEnd = _vecd.end();
    for (; dIt != dEnd; ++dIt) {
        totalVariance += pow((*dIt - mean),2);
    }
    double variance = totalVariance / _vecd.size();
    return sqrt(variance);
}

double VectorOps::GetSstddev(void) 
{
    double mean = GetMean();
    // get the variance
    double totalVariance = 0.0;
    vector<double>::const_iterator dIt  = _vecd.begin();
    vector<double>::const_iterator dEnd = _vecd.end();
    for (; dIt != dEnd; ++dIt) {
        totalVariance += pow((*dIt - mean),2);
    }
    double variance = totalVariance / (_vecd.size() - 1);
    return sqrt(variance);
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


string VectorOps::GetMode(void)
{
    string mode;
    int max_count = 0;
    int curr_count = 0;
    
    sort(_vecs.begin(), _vecs.end());

    string prev = "";
    vector<string>::const_iterator it = _vecs.begin();
    while (it != _vecs.end())
    {
        if (*it != prev && prev != "")
        {
            if (curr_count > max_count)
            {
                max_count = curr_count;
                mode = prev;
            }
            curr_count = 1;
        }
        else { curr_count++; }
        prev = *it;
        ++it;
    }
    if (curr_count > max_count)
        mode = prev;
    return mode;
}


string VectorOps::GetAntiMode(void)
{
    string antimode;
    int min_count = INT_MAX;
    int curr_count = 0;
    
    sort(_vecs.begin(), _vecs.end());

    string prev = "";
    vector<string>::const_iterator it = _vecs.begin();
    while (it != _vecs.end())
    {
        if (*it != prev && prev != "")
        {
            if (curr_count < min_count)
            {
                min_count = curr_count;
                antimode = prev;
            }
            curr_count = 1;
        }
        else { curr_count++; }
        prev = *it;
        ++it;
    }
    if (curr_count < min_count)
        antimode = prev;
    return antimode;
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

string VectorOps::GetConcat(void)
{
    ostringstream collapse;
    for( size_t i = 0; i < _vecs.size(); i++ )
        collapse << _vecs[i];
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

string VectorOps::GetFreqDesc(void)
{
    // compute the frequency of each unique value
    map<string, int> freqs;
    vector<string>::const_iterator dIt  = _vecs.begin();
    vector<string>::const_iterator dEnd = _vecs.end();
    for (; dIt != dEnd; ++dIt) {
        freqs[*dIt]++;
    }
    // pair for the num times a values was
    // observed (1) and the value itself (2)
    pair<int, string> freqPair;
    vector< pair<int, string> > freqList;

    // create a list of pairs of all the observed values (second)
    // and their occurences (first)
    map<string,int>::const_iterator mapIter = freqs.begin();
    map<string,int>::const_iterator mapEnd  = freqs.end();
    for(; mapIter != mapEnd; ++mapIter)
        freqList.push_back( make_pair(mapIter->second, mapIter->first) );

    // sort the list of pairs in the requested order by the frequency
    // this will make the value that was observed least/most bubble to the top
    sort(freqList.begin(), freqList.end(), ValueGreaterThan());

    // record all of the values and their frequencies.
    ostringstream buffer;
    vector< pair<int, string> >::const_iterator iter    = freqList.begin();
    vector< pair<int, string> >::const_iterator iterEnd = freqList.end();
    for (; iter != iterEnd; ++iter)
        buffer << iter->second << ":" << iter->first << ",";
    
    return buffer.str();
}


string VectorOps::GetFreqAsc(void)
{
    // compute the frequency of each unique value
    map<string, int> freqs;
    vector<string>::const_iterator dIt  = _vecs.begin();
    vector<string>::const_iterator dEnd = _vecs.end();
    for (; dIt != dEnd; ++dIt) {
        freqs[*dIt]++;
    }
    // pair for the num times a values was
    // observed (1) and the value itself (2)
    pair<int, string> freqPair;
    vector< pair<int, string> > freqList;

    // create a list of pairs of all the observed values (second)
    // and their occurences (first)
    map<string,int>::const_iterator mapIter = freqs.begin();
    map<string,int>::const_iterator mapEnd  = freqs.end();
    for(; mapIter != mapEnd; ++mapIter)
        freqList.push_back( make_pair(mapIter->second, mapIter->first) );

    // sort the list of pairs in the requested order by the frequency
    // this will make the value that was observed least/most bubble to the top
    sort(freqList.begin(), freqList.end(), ValueLessThan());

    // record all of the values and their frequencies.
    ostringstream buffer;
    vector< pair<int, string> >::const_iterator iter    = freqList.begin();
    vector< pair<int, string> >::const_iterator iterEnd = freqList.end();
    for (; iter != iterEnd; ++iter)
        buffer << iter->second << ":" << iter->first << ",";
    
    return buffer.str();
}
