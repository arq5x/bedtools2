/*****************************************************************************
groupBy.cpp

(c) 2009, 2010, 2011 - Aaron Quinlan
Center for Public Health Genomics
University of Virginia
aaronquinlan@gmail.com

Licenced under the MIT license.
******************************************************************************/
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <exception>
#include <stdexcept> // out_of_range exception

#include "version.h"
#include "lineFileUtilities.h"
#include "tabFile.h"
using namespace std;


const int PRECISION = 21;


// define our program name
#define PROGRAM_NAME "bedtools groupby"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))
#define LOOKS_LIKE_A_PARAM(string) (strlen(string)>0 && string[0]=='-')


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

// function declarations
void groupby_help(void);
void GroupBy(const string &inFile, const vector<int> &groupColumns, const vector<int> &opColumns, const vector<string> &ops, const bool printOriginalLine, const bool printHeaderLine, const bool InputHaveHeaderLine, const bool ignoreCase);
void PrintHeaderLine(const vector<string> &InputFields, const vector<int> &groupColumns, const vector<int> &opColumns, const vector<string> &ops, const bool PrintFullInputLine, const bool InputHaveHeaderLine);
void ReportSummary(const vector<string> &group, const vector<vector<string> > &data, const vector<string> &ops);
void addValue (const vector<string> &fromList, vector<string> &toList, int index, int lineNum, const bool ignoreCase);
float ToFloat (string element);
double ToDouble(const string &element);
void TabPrintPost (string element);
void TabPrintPre (string element);
void CommaPrint (string element);

int groupby_main(int argc, char* argv[]) {

    // input files
    string inFile             = "stdin";
    string groupColumnsString = "1,2,3";
    string opsColumnString;
    string opsString;

    // our configuration variables
    bool showHelp          = false;
    bool haveOpColumns     = false;
    bool haveOps           = true;
    bool printOriginalLine = false;
    bool printHeaderLine   = false;
    bool InputHaveHeaderLine = false;
    bool ignoreCase    = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) groupby_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) >= argc || LOOKS_LIKE_A_PARAM(argv[i+1])) {
                cerr << endl << "*****ERROR: -i parameter requires a value." << endl << endl;
                groupby_help();
                break;
            }
            else {
                inFile     = argv[i + 1];
                i++;
            }
        }
        else if (PARAMETER_CHECK("-grp", 4, parameterLength) || PARAMETER_CHECK("-g", 2, parameterLength)) {
            if ((i+1) >= argc || LOOKS_LIKE_A_PARAM(argv[i+1])) {
                cerr << endl << "*****ERROR: -grp parameter requires a value." << endl << endl;
                groupby_help();
                break;
            }
            else {
                groupColumnsString     = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-opCols", 7, parameterLength) || PARAMETER_CHECK("-c", 2, parameterLength)) {
            if ((i+1) >= argc || LOOKS_LIKE_A_PARAM(argv[i+1])) {
                cerr << endl << "*****ERROR: -opCols parameter requires a value." << endl << endl;
                groupby_help();
                break;
            }
            else {
                haveOpColumns       = true;
                opsColumnString     = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-ops", 4, parameterLength) || PARAMETER_CHECK("-o", 2, parameterLength)) {
            if ((i+1) >= argc || LOOKS_LIKE_A_PARAM(argv[i+1])) {
                cerr << endl << "*****ERROR: -ops parameter requires a value." << endl << endl;
                groupby_help();
                break;
            }
            else {
                haveOps    = true;
                opsString  = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-full", 5, parameterLength)) {
            printOriginalLine = true;
        }
        else if(PARAMETER_CHECK("-outheader", 10, parameterLength)) {
            printHeaderLine = true;
        }
        else if(PARAMETER_CHECK("-inheader", 9, parameterLength)) {
            InputHaveHeaderLine = true;
        }
        else if(PARAMETER_CHECK("-header", 7, parameterLength)) {
            InputHaveHeaderLine = true;
            printHeaderLine = true;
        }
        else if(PARAMETER_CHECK("-ignorecase", 11, parameterLength)) {
            ignoreCase = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!haveOpColumns) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -opCols." << endl << "*****" << endl;
        showHelp = true;
    }
    // split the opsString into discrete operations and make sure they are all valid.
    vector<string> ops;
    opsString.erase(remove_if(opsString.begin(),opsString.end(),::isspace),opsString.end());
    Tokenize(opsString, ops, ",");
    for( size_t i = 0; i < ops.size(); i++ ) {
        if ((ops[i] != "sum")  && (ops[i] != "max")    && (ops[i] != "min") && (ops[i] != "mean") &&
            (ops[i] != "mode") && (ops[i] != "median") && (ops[i] != "antimode") && (ops[i] != "stdev") &&
            (ops[i] != "sstdev") && (ops[i] != "count") && (ops[i] != "collapse") && (ops[i] != "distinct") &&
            (ops[i] != "concat") && (ops[i] != "freqdesc") && (ops[i] != "freqasc")) 
        {
            cerr << endl << "*****" << endl << "*****ERROR: Invalid operation selection \"" << ops[i] << endl << "\"  *****" << endl;
            showHelp = true;
        }
    }
    if (!showHelp) {

        // Split the column string sent by the user into discrete column numbers
        // A comma separated string is expected.
        vector<int> groupColumnsInt;
        Tokenize(groupColumnsString, groupColumnsInt, ",");

        vector<int> opColumnsInt;
        Tokenize(opsColumnString, opColumnsInt, ",");

        // sanity check the group columns
        for(size_t i = 0; i < groupColumnsInt.size(); ++i) {
            int groupColumnInt = groupColumnsInt[i];
            if (groupColumnInt < 1) {
                cerr << endl << "*****" << endl << "*****ERROR: group columns must be >=1. " << endl << "*****" << endl;
                groupby_help();
            }
        }

        // sanity check the op columns
        for(size_t i = 0; i < opColumnsInt.size(); ++i) {
            int opColumnInt = opColumnsInt[i];
            if (opColumnInt < 1) {
                cerr << endl << "*****" << endl << "*****ERROR: op columns must be >=1. " << endl << "*****" << endl;
                groupby_help();
            }
        }

        // sanity check that there are equal number of opColumns and ops
        if (ops.size() != opColumnsInt.size()) {
            cerr << endl << "*****" << endl << "*****ERROR: There must be equal number of ops and opCols. " << endl << "*****" << endl;
            groupby_help();
        }
        GroupBy(inFile, groupColumnsInt, opColumnsInt, ops,
            printOriginalLine, printHeaderLine, InputHaveHeaderLine,
            ignoreCase);
    }
    else {
        groupby_help();
    }
    return 0;
}

void groupby_help(void) {

    cerr << "\nTool:    bedtools intersect (aka intersectBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Summarizes a dataset column based upon" << endl;
    cerr << "\t common column groupings. Akin to the SQL \"group by\" command." << endl << endl;

    cerr << "Usage:\t " << PROGRAM_NAME << " -g [group_column(s)] -c [op_column(s)] -o [ops] " << endl;
    cerr << "\t "     << "cat [FILE] | " << PROGRAM_NAME << " -g [group_column(s)] -c [op_column(s)] -o [ops] " << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-i\t\t"        << "Input file. Assumes \"stdin\" if omitted." << endl << endl;

    cerr << "\t-g -grp\t\t"      << "Specify the columns (1-based) for the grouping." << endl;
    cerr                         << "\t\t\tThe columns must be comma separated." << endl;
    cerr                         << "\t\t\t- Default: 1,2,3" << endl << endl;

    cerr << "\t-c -opCols\t"   << "Specify the column (1-based) that should be summarized." << endl;
    cerr                         << "\t\t\t- Required." << endl << endl;

    cerr << "\t-o -ops\t\t"      << "Specify the operation that should be applied to opCol." << endl;
    cerr                         << "\t\t\tValid operations:" << endl;
    cerr                         << "\t\t\t    sum, count, min, max," << endl;
    cerr                         << "\t\t\t    mean, median, mode, antimode," << endl;
    cerr                         << "\t\t\t    stdev, sstdev (sample standard dev.)," << endl;
    cerr                         << "\t\t\t    collapse (i.e., print a comma separated list (duplicates allowed)), " << endl;
    cerr                         << "\t\t\t    distinct (i.e., print a comma separated list (NO duplicates allowed)), " << endl;
    cerr                         << "\t\t\t    concat   (i.e., merge values into a single, non-delimited string), " << endl;
    cerr                         << "\t\t\t    freqdesc (i.e., print desc. list of values:freq)" << endl;
    cerr                         << "\t\t\t    freqasc (i.e., print asc. list of values:freq)" << endl;
    cerr                         << "\t\t\t- Default: sum" << endl << endl;

    cerr << "\t-full\t\t"   << "Print all columns from input file." << endl;
    cerr            << "\t\t\tDefault: print only grouped columns." << endl << endl;

    cerr << "\t-inheader\t" << "Input file has a header line - the first line will be ignored." << endl << endl ;

    cerr << "\t-outheader\t"    << "Print header line in the output, detailing the column names. " << endl;
    cerr            << "\t\t\tIf the input file has headers (-inheader), the output file" << endl;
    cerr            << "\t\t\twill use the input's column names." << endl;
    cerr            << "\t\t\tIf the input file has no headers, the output file" << endl;
    cerr            << "\t\t\twill use \"col_1\", \"col_2\", etc. as the column names." << endl << endl;

    cerr << "\t-header\t\t" << "same as '-inheader -outheader'" << endl << endl;

    cerr << "\t-ignorecase\t"   << "Group values regardless of upper/lower case." << endl << endl;

    cerr << "Examples: " << endl;
    cerr << "\t$ cat ex1.out" << endl;
    cerr << "\tchr1 10  20  A   chr1    15  25  B.1 1000    ATAT" << endl;
    cerr << "\tchr1 10  20  A   chr1    25  35  B.2 10000   CGCG" << endl << endl;
    cerr << "\t$ groupBy -i ex1.out -g 1,2,3,4 -c 9 -o sum" << endl;
    cerr << "\tchr1 10  20  A   11000" << endl << endl;
    cerr << "\t$ groupBy -i ex1.out -grp 1,2,3,4 -opCols 9,9 -ops sum,max" << endl;
    cerr << "\tchr1 10  20  A   11000   10000" << endl << endl;
    cerr << "\t$ groupBy -i ex1.out -g 1,2,3,4 -c 8,9 -o collapse,mean" << endl;
    cerr << "\tchr1 10  20  A   B.1,B.2,    5500" << endl << endl;
    cerr << "\t$ cat ex1.out | groupBy -g 1,2,3,4 -c 8,9 -o collapse,mean" << endl;
    cerr << "\tchr1 10  20  A   B.1,B.2,    5500" << endl << endl;
    cerr << "\t$ cat ex1.out | groupBy -g 1,2,3,4 -c 10 -o concat" << endl;
    cerr << "\tchr1 10  20  A   ATATCGCG" << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  The input file/stream should be sorted/grouped by the -grp. columns" << endl;
    cerr << "\t(2)  If -i is unspecified, input is assumed to come from stdin." << endl << endl;


    // end the program here
    exit(1);

}


void GroupBy (const string &inFile,
    const vector<int> &groupColumns,
    const vector<int> &opColumns,
    const vector<string> &ops,
    const bool printOriginalLine,
    const bool printHeaderLine,
    const bool InputHaveHeaderLine,
    const bool ignoreCase) {

    // current line number
    int lineNum = 0;
    // string representing current line
    string inLine;

    // vector of strings holding the tokenized current line
    vector<string>  inFields;
    vector<string>  inFieldsFirstLineInGroup;
    inFields.reserve(20);

    // keys for the current and previous group
    vector<string>  prevGroup(0);
    vector<string>  currGroup(0);

    // vector (one per column) of vector (one per value/column) of the opColumn values for the current group
    vector< vector<string> >  values;
    for( size_t i = 0; i < opColumns.size(); i++ ) {
        values.push_back( vector<string>() );
    }

    bool    first_line = true;

    // check the status of the current line
    TabLineStatus tabLineStatus;

    // open a new tab file, loop through it line by line
    // and summarize the data for a given group when the group
    // fields change
    TabFile *_tab = new TabFile(inFile);
    _tab->Open();
    while ((tabLineStatus = _tab->GetNextTabLine(inFields, lineNum)) != TAB_INVALID) {
        if (tabLineStatus == TAB_VALID) {

            if (first_line) {
                first_line = false;
                if (printHeaderLine)
                    PrintHeaderLine(inFields, groupColumns, opColumns, ops,
                    printOriginalLine, InputHaveHeaderLine);

                if (InputHaveHeaderLine) {
                    inFields.clear();
                    continue; // no need to process this line - it's the header
                }
            }

            if (inFieldsFirstLineInGroup.empty()) //first line in file? - save it
                inFieldsFirstLineInGroup = inFields;

            // build the group vector for the current line
            currGroup.clear();
            vector<int>::const_iterator gIt  = groupColumns.begin();
            vector<int>::const_iterator gEnd = groupColumns.end();
            for (; gIt != gEnd; ++gIt)
                addValue(inFields, currGroup, (*gIt-1), lineNum, ignoreCase);

            // there has been a group change
            if ((currGroup != prevGroup) && (prevGroup.size() > 0)) {
                // Summarize this group
                ReportSummary(printOriginalLine?inFieldsFirstLineInGroup:prevGroup, values, ops);
                // reset and add the first value for the next group.
                values.clear();
                for( size_t i = 0; i < opColumns.size(); i++ ) {
                    values.push_back( vector<string>() );
                    addValue(inFields, values[i], opColumns[i]-1, lineNum, ignoreCase);
                }
                inFieldsFirstLineInGroup = inFields;
            }
            // we're still dealing with the same group
            else {
                for( size_t i = 0; i < opColumns.size(); i++ )
                    addValue(inFields, values[i], opColumns[i]-1, lineNum, ignoreCase);
            }
            // reset for the next line
            prevGroup = currGroup;
        }
        inFields.clear();
    }
    // report the last group
    ReportSummary(printOriginalLine?inFieldsFirstLineInGroup:currGroup, values, ops);
    _tab->Close();
}


void ReportSummary(const vector<string> &group, const vector<vector<string> > &data, const vector<string> &ops) {

    vector<string> result;
    for( size_t i = 0; i < data.size(); i++ ) {

        string op = ops[i];
        std::stringstream buffer;
        vector<double> dataF;
        // are we doing a numeric conversion?  if so, convert the strings to doubles.
        if ((op == "sum") || (op == "max") || (op == "min") || (op == "mean") ||
            (op == "median") || (op == "stdev") || (op == "sstdev"))
        {
            transform(data[i].begin(), data[i].end(), back_inserter(dataF), ToDouble);
        }

        if (op == "sum") {
            // sum them up
            double total = accumulate(dataF.begin(), dataF.end(), 0.0);
            buffer << setprecision (PRECISION) << total;
            result.push_back(buffer.str());
        }
        else if (op == "collapse") {
            string collapse;
            for( size_t j = 0; j < data[i].size(); j++ ) {//Ugly, but cannot use back_inserter
                if (j>0)
                    collapse.append(",");
                collapse.append(data[i][j]);
            }
            result.push_back(collapse);
        }
        else if (op == "distinct") {
            string distinct;
            // get the current column's data
            vector<string> col_data = data[i];
            // remove duplicate entries from the vector
            // http://stackoverflow.com/questions/1041620/most-efficient-way-to-erase-duplicates-and-sort-a-c-vector
            sort( col_data.begin(), col_data.end() );
            col_data.erase( unique( col_data.begin(), col_data.end() ), col_data.end() );
            
            for( size_t j = 0; j < col_data.size(); j++ ) {//Ugly, but cannot use back_inserter
                if (j>0)
                    distinct.append(",");
                distinct.append(col_data[j]);
            }
            result.push_back(distinct);
        }
        else if (op == "concat") {
            string concat;
            for( size_t j = 0; j < data[i].size(); j++ ) {//Ugly, but cannot use back_inserter
                concat.append(data[i][j]);
            }
            result.push_back(concat);
        }
        else if (op == "min") {
            buffer << setprecision (PRECISION) << *min_element( dataF.begin(), dataF.end() );
            result.push_back(buffer.str());
        }
        else if (op == "max") {
            buffer << setprecision (PRECISION) << *max_element( dataF.begin(), dataF.end() );
            result.push_back(buffer.str());
        }
        else if (op == "mean") {
            double total = accumulate(dataF.begin(), dataF.end(), 0.0);
            double mean = total / dataF.size();
            buffer << setprecision (PRECISION) << mean;
            result.push_back(buffer.str());
        }
        else if (op == "median") {
            double median = 0.0;
            sort(dataF.begin(), dataF.end());
            int totalLines = dataF.size();
            if ((totalLines % 2) > 0) {
                long mid;
                mid = totalLines / 2;
                median = dataF[mid];
            }
            else {
                long midLow, midHigh;
                midLow = (totalLines / 2) - 1;
                midHigh = (totalLines / 2);
                median = (dataF[midLow] + dataF[midHigh]) / 2.0;
            }
            buffer << setprecision (PRECISION) << median;
            result.push_back(buffer.str());
        }
        else if (op == "count") {
            buffer << setprecision (PRECISION) << data[i].size();
            result.push_back(buffer.str());
        }
        else if ((op == "mode") || (op == "antimode") ||
        (op == "freqdesc") || (op == "freqasc")) {
            // compute the frequency of each unique value
            map<string, int> freqs;
            vector<string>::const_iterator dIt  = data[i].begin();
            vector<string>::const_iterator dEnd = data[i].end();
            for (; dIt != dEnd; ++dIt) {
                freqs[*dIt]++;
            }

            // grab the mode and the anti mode
            string mode, antiMode;
            int    count = 0;
            int minCount = INT_MAX;
            for(map<string,int>::const_iterator iter = freqs.begin(); iter != freqs.end(); ++iter) {
                if (iter->second > count) {
                    mode = iter->first;
                    count = iter->second;
                }
                if (iter->second < minCount) {
                    antiMode = iter->first;
                    minCount = iter->second;
                }
            }
            // report
            if (op == "mode") {
                buffer << setprecision (PRECISION) << mode;
                result.push_back(buffer.str());
            }
            else if (op == "antimode") {
                buffer << setprecision (PRECISION) << antiMode;
                result.push_back(buffer.str());
            }
            else if (op == "freqdesc" || op == "freqasc") {
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
                if (op == "freqdesc")
                    sort(freqList.begin(), freqList.end(), ValueGreaterThan());
                else if (op == "freqasc")
                    sort(freqList.begin(), freqList.end(), ValueLessThan());

                // record all of the values and their frequencies.
                vector< pair<int, string> >::const_iterator iter    = freqList.begin();
                vector< pair<int, string> >::const_iterator iterEnd = freqList.end();
                for (; iter != iterEnd; ++iter)
                    buffer << iter->second << ":" << iter->first << ",";
                result.push_back(buffer.str());
            }
        }
        else if (op == "stdev" || op == "sstdev") {
            // get the mean
            double total = accumulate(dataF.begin(), dataF.end(), 0.0);
            double mean = total / dataF.size();
            // get the variance
            double totalVariance = 0.0;
            vector<double>::const_iterator dIt  = dataF.begin();
            vector<double>::const_iterator dEnd = dataF.end();
            for (; dIt != dEnd; ++dIt) {
                totalVariance += pow((*dIt - mean),2);
            }
            double variance = 0.0;
            if (op == "stdev") {
                variance = totalVariance / dataF.size();
            }
            else if (op == "sstdev" && dataF.size() > 1) {
                variance = totalVariance / (dataF.size() - 1);
            }
            double stddev = sqrt(variance);
            // report
            buffer << setprecision (PRECISION) << stddev;
            result.push_back(buffer.str());
        }
    }
    for_each(group.begin(), group.end(), TabPrintPost);
    cout << *result.begin();
    for_each(++result.begin(), result.end(), TabPrintPre);
    cout << endl; //Gets rid of extraneous tab
}


void addValue (const vector<string> &fromList, vector<string> &toList, int index, int lineNum, const bool ignoreCase) {
    try {
        string s(fromList.at(index));
        if(ignoreCase)
            transform(s.begin(),s.end(),s.begin(),::tolower);
        toList.push_back(s);
    }
    catch(std::out_of_range& e) {
        cerr << endl << "*****" << endl << "*****ERROR: requested column exceeds the number of columns in file at line "
            << lineNum << ". Exiting." << endl << "*****" << endl;
        exit(1);
    }
}


float ToFloat (string element) {
    return atof(element.c_str());
}

void TabPrintPost (string element) {
    cout << element << "\t";
}

void TabPrintPre (string element) {
    cout << "\t" << element;
}

void CommaPrint (string element) {
    cout << element << ",";
}

double ToDouble(const string &element) {
    std::istringstream i(element);
    double x;
    if (!(i >> x)) {
        cerr << "Error: Could not properly convert string to numeric (\"" + element + "\")" << endl;
        exit(1);
    }
    return x;
}

inline string ColumnHeaderName(const vector<string> &inFields, const size_t FieldIndex,
    bool InputHaveHeaderLine)
{
    stringstream s;
    if (InputHaveHeaderLine)
        s << inFields[FieldIndex-1];
    else
        s << "col_" << (FieldIndex);
    return s.str();
}

void PrintHeaderLine(const vector<string> &inFields,
    const vector<int> &groupColumns,
    const vector<int> &opColumns,
    const vector<string> &ops,
    const bool PrintFullInputLine,
    const bool InputHaveHeaderLine)
{
    vector<string> header;

    //Header fields of input file
    if (PrintFullInputLine) {
        //All input columns
        for (size_t i=0;i<inFields.size();++i)
            header.push_back( ColumnHeaderName(inFields, i+1, InputHaveHeaderLine) );
    } else {
        //Only the columns that are actually used in the grouped operations
        for (size_t i=0;i<groupColumns.size();++i)
            header.push_back( ColumnHeaderName(inFields, groupColumns[i], InputHaveHeaderLine) );
    }

    //Header fields of output columns, by operation
    for (size_t i=0; i<opColumns.size();++i) {
        stringstream s;
        s << ops[i] << "(" << ColumnHeaderName(inFields, opColumns[i], InputHaveHeaderLine) << ")";
        header.push_back(s.str());
    }

    //print Header Line
    for (size_t i=0; i<header.size();++i)
        cout << header[i] << ((i<header.size()-1)?"\t":"\n");
}
