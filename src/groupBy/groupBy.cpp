/*****************************************************************************
  overlap.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include <vector>
#include <map>
#include <numeric>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include "version.h"
#include "lineFileUtilities.h"
#include "tabFile.h"
using namespace std;


// define our program name
#define PROGRAM_NAME "groupBy"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);
void GroupBy(const string &inFile, const vector<int> &groupColumns, int opColumn, string op);
void ReportSummary(const vector<string> &group, const vector<string> &data, string op);
float ToFloat (string element);
void TabPrint (string element);
void CommaPrint (string element);
            
int main(int argc, char* argv[]) {

	// input files
	string inFile;
	string groupColumnsString = "1,2,3";
    string opColumnString;
    string op = "sum";
    
	// our configuration variables
	bool showHelp         = false;
	bool haveInFile       = false;
	bool haveGroupColumns = false;
	bool haveOpColumn     = false;	
	bool haveOp           = true;	
	
	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
		(PARAMETER_CHECK("--help", 5, parameterLength))) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();

	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-i", 2, parameterLength)) {
			if ((i+1) < argc) {
				haveInFile = true;
				inFile     = argv[i + 1];
				i++;
			}
		}
		else if(PARAMETER_CHECK("-grp", 4, parameterLength)) {
			haveGroupColumns       = true;
			groupColumnsString     = argv[i + 1];
			i++;			
		}
		else if(PARAMETER_CHECK("-opCol", 6, parameterLength)) {
			haveOpColumn       = true;
			opColumnString     = argv[i + 1];
			i++;			
		}
		else if(PARAMETER_CHECK("-op", 3, parameterLength)) {
			haveOp = true;
			op     = argv[i + 1];
			i++;			
		}
		else {
			cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have an input files
	if (!haveInFile) {
		cerr << endl << "*****" << endl << "*****ERROR: Need -i file. " << endl << "*****" << endl;
		showHelp = true;
	}
	if (!haveOpColumn) {
		cerr << endl << "*****" << endl << "*****ERROR: Need -opCol." << endl << "*****" << endl;
		showHelp = true;
	}
	if ((op != "sum")  && (op != "max")    && (op != "min") && (op != "mean") &&
	    (op != "mode") && (op != "median") && (op != "antimode") && (op != "stdev") &&
	    (op != "count") && (op != "collapse")) {
		    cerr << endl << "*****" << endl << "*****ERROR: Invalid op selection. " << endl << "*****" << endl;	
		    showHelp = true;		                    
	}
	if (!showHelp) {
		
		// Split the column string sent by the user into discrete column numbers
		// A comma separated string is expected.
		vector<int> groupColumnsInt;
		Tokenize(groupColumnsString, groupColumnsInt, ",");
		
        int opColumnInt = atoi(opColumnString.c_str());

		//DetermineInput(inFile, groupColumnsInt, opColumnInt, op);
		GroupBy(inFile, groupColumnsInt, opColumnInt, op);
	}	
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {

	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Summarizes a dataset column based upon" << endl;
	cerr << "\t common column groupings. Akin to the SQL \"group by\" command." << endl << endl;
	
	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <input> -opCol <> " << endl << endl;

	cerr << "Options: " << endl;
	cerr << "\t-i\t"	    << "Input file. Use \"stdin\" for pipes." << endl << endl;
	
	cerr << "\t-grp\t"	    << "Specify the columns (1-based) for the grouping." << endl;
	cerr 					<< "\t\tThe columns must be comma separated." << endl;
	cerr					<< "\t\tDefault: 1,2,3" << endl << endl;	

	cerr << "\t-opCol\t"	<< "Specify the column (1-based) that should be summarized." << endl;
	cerr 					<< "\t\tRequired." << endl << endl;

	cerr << "\t-op\t"	    << "Specify the operation that should be applied to opCol." << endl;
	cerr 					<< "\t\tValid operations: sum, count, min, max, mean, median," << endl;
    cerr                    << "\t\tmode, antimode, stdev, collapse (i.e., print a comma separated list)" << endl;
    cerr                    << "\t\tDefault: sum" << endl;

	cerr << "Examples: " << endl;
	cerr << "\t$ cat test.out" << endl;
	cerr << "\tchr1	10	20	A	chr1	15	25	1000" << endl;
	cerr << "\tchr1	10	20	A	chr1	25	35	10000" << endl << endl;
	cerr << "\t$ cat test.out | groupBy -i stdin -grpCols 1,2,3,4 -opCol 8 -op sum" << endl;
	cerr << "\tchr1	10	20	A	11000" << endl << endl;
	cerr << "\t$ cat test.out | groupBy -i stdin -grpCols 1,2,3,4 -opCol 8 -op max" << endl;
	cerr << "\tchr1	10	20	A	1000" << endl << endl;
	cerr << "\t$ cat test.out | groupBy -i stdin -grpCols 1,2,3,4 -opCol 8 -op mean" << endl;
	cerr << "\tchr1	10	20	A	5500" << endl << endl;
	
	cerr << "Notes: " << endl;
	cerr << "\t(1)  The input file/stream should be sorted/grouped by the -grpCols." << endl << endl;

	
	// end the program here
	exit(1);

}


void GroupBy(const string &inFile, const vector<int> &groupColumns, int opColumn, string op) {
	
	// current line number
	int lineNum = 0;
	// string representing current line
	string inLine;
	
	// vector of strings holding the tokenized current line
	vector<string>  inFields;
    inFields.reserve(20);
    
    // keys for the current and previous group
    vector<string>  prevGroup(0);
    vector<string>  currGroup(0);
    
    // vector of the opColumn values for the current group
    vector<string>  values;
    values.reserve(100000);
    
    // check the status of the current line
    TabLineStatus tabLineStatus;
    
    // open a new tab file, loop through it line by line
    // and summarize the data for a given group when the group
    // fields change
    TabFile *_tab = new TabFile(inFile);    
    _tab->Open();
	while ((tabLineStatus = _tab->GetNextTabLine(inFields, lineNum)) != TAB_INVALID) {
		if (tabLineStatus == TAB_VALID) {
		    // build the group vector for the current line
		    currGroup.clear();
            vector<int>::const_iterator gIt  = groupColumns.begin();
            vector<int>::const_iterator gEnd = groupColumns.end();
            for (; gIt != gEnd; ++gIt) {currGroup.push_back(inFields[*gIt-1]);}

            // group change
            if ((currGroup != prevGroup) && (prevGroup.size() > 0)) {
                // Summarize this group
                ReportSummary(prevGroup, values, op);
                values.clear();
                values.push_back(inFields[opColumn-1].c_str());
            }
            // same group
            else
                values.push_back(inFields[opColumn-1].c_str());

            // reset for the next line
            prevGroup = currGroup;
            inFields.clear();
	    }
    }
    // report the last group
    values.clear();
    values.push_back(inFields[opColumn-1].c_str());
    ReportSummary(currGroup, values, op);
    
    _tab->Close();
}


void ReportSummary(const vector<string> &group, const vector<string> &data, string op) {
    vector<double> dataF;
    // convert to floats
    transform(data.begin(), data.end(), back_inserter(dataF), ToFloat);
    
    if (op == "sum") {        
        // sum them up
        double total = accumulate(dataF.begin(), dataF.end(), 0.0);
        for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << total << endl; 
    }
    if (op == "collapse") {        
        for_each(group.begin(), group.end(), TabPrint);
        for_each(data.begin(), data.end(), CommaPrint);
        cout << endl;
    }
    else if (op == "min") {
        for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << *min_element( dataF.begin(), dataF.end() ) << endl;
    }
    else if (op == "max") {
        for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << *max_element( dataF.begin(), dataF.end() ) << endl;
    }
    else if (op == "mean") {
        double total = accumulate(dataF.begin(), dataF.end(), 0.0);
        double mean = total / dataF.size();
        for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << mean << endl;
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
		for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << median << endl;
    }
    else if (op == "count") {
        for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << data.size() << endl;
    }
    else if ((op == "mode") || (op == "antimode")) {
        // compute the frequency of each unique value
        map<string, int> freqs;
        vector<string>::const_iterator dIt  = data.begin();
        vector<string>::const_iterator dEnd = data.end();
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
        for_each(group.begin(), group.end(), TabPrint);
        if (op == "mode")
            cout << setprecision (7) << mode << endl;
        else if (op == "antimode")
            cout << setprecision (7) << antiMode << endl;
    }
    else if (op == "stdev") {
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
        double variance = totalVariance / dataF.size();
        double stddev = sqrt(variance);
        // report
        for_each(group.begin(), group.end(), TabPrint);
        cout << setprecision (7) << stddev << endl;
    }
}


float ToFloat (string element) {
    return atof(element.c_str());
}

void TabPrint (string element) {
    cout << element << "\t";
}

void CommaPrint (string element) {
    cout << element << ",";
}
