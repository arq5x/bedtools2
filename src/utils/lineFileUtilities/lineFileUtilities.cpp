// 
//  lineFileUtilities.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Contains common functions for processing text files.
//

#include "lineFileUtilities.h"
//***********************************************
// lineFileUtilities:
// 		Common Functions
//***********************************************

void Tokenize(const string& str, vector<string>& tokens)
{

/*
	string::size_type start = 0;
	string::size_type end = 0;

	cout << str << endl; 
	while ((start = str.find_first_not_of(" ", start)) != string::npos) {
	    end = str.find_first_of(" ", start);
		cout << start << ":" << end << endl;
	    tokens.push_back(str.substr(start, end - start));
	    start = end;
	  }
*/

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of("\t", 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of("\t", lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of("\t", pos);
		// Find next "non-delimiter"
		pos = str.find_first_of("\t", lastPos);
	}

}


