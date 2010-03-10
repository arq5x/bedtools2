// 
//  lineFileUtilities.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Contains common functions for processing text files.
//
#include <sstream>
#include <iostream>
#include "lineFileUtilities.h"
//***********************************************
// lineFileUtilities:
// 		Common Functions
//***********************************************
void Tokenize(string str, vector<string> &tokens, const string &delimiter) {

	/*
	//method to tokenize on any whitespace
    string buf; 			// Have a buffer string
    stringstream ss(str); // Insert the string into a stream
    while (ss >> buf)
		//cout << buf << endl;
        tokens.push_back(buf);
	*/

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiter, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter, lastPos);
	}
	
}


