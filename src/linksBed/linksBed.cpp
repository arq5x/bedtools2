// 
//  linksBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Creates HTML or UCSC browser coordinates from a BED file.
//
#include "lineFileUtilities.h"
#include "linksBed.h"

//
// Constructor
//
BedLinks::BedLinks(string &bedFile, string &base, string &org, string &db) {
	this->bedFile = bedFile;
	this->bed = new BedFile(bedFile);
	
	this->base = base;
	this->org = org;
	this->db = db;
}

//
// Destructor
//
BedLinks::~BedLinks(void) {
}


void BedLinks::WriteURL(BED &bed, string &base) {

	string position = bed.chrom;
	std::stringstream posStream;
	posStream << ":" << bed.start << "-" << bed.end;
	position.append(posStream.str());

	cout << "<tr>" << endl;
		cout << "\t<td>" << endl;
			cout << "\t\t<a href=" << base << position << ">";
			cout << bed.chrom << ":" << bed.start << "-" << bed.end; 
			cout << "</a>" << endl;
		cout << "\t</td>" << endl;	

		if (this->bed->bedType == 4) {
			cout << "\t<td>" << endl;
			cout << bed.name << endl;
			cout << "\t</td>" << endl;
		}
		else if (this->bed->bedType == 5) {
			cout << "\t<td>" << endl;
			cout << bed.name << endl;
			cout << "\t</td>" << endl;
			
			cout << "\t<td>" << endl;
			cout << bed.score << endl;
			cout << "\t</td>" << endl;
		}
		else if (this->bed->bedType == 6) {
			cout << "\t<td>" << endl;
			cout << bed.name << endl;
			cout << "\t</td>" << endl;
			
			cout << "\t<td>" << endl;
			cout << bed.score << endl;
			cout << "\t</td>" << endl;
			
			cout << "\t<td>" << endl;
			cout << bed.strand << endl;
			cout << "\t</td>" << endl;
		}		
		cout << "</tr>" << endl;
}


void BedLinks::LinksBed() {


	// construct the html base.
	string org = this->org;
	string db = this->db;
	string base = this->base;
	base.append("/cgi-bin/hgTracks?org=");
	base.append(org);
	base.append("&db=");
	base.append(db);
	base.append("&position=");
	
	// create the HTML header
	cout << "<html>" << endl <<"\t<body>" << endl; 
	cout << "<title>" << this->bedFile << "</title>" << endl;
	
	// start the table of entries
	cout << "<br>Firefox users: Press and hold the \"apple\" or \"alt\" key and click link to open in new tab." << endl;
	cout << "<p style=\"font-family:courier\">" << endl;
	cout << "<table border=\"0\" align=\"justify\"" << endl;
	
	
	// Are we dealing with a BED file or a BED passed via stdin?

	// Case 1: Proper BED File.
	if ( (this->bedFile != "") && (this->bedFile != "stdin") ) {

		// open the BED file for reading                                                                                                                                      
		ifstream bed(bedFile.c_str(), ios::in);
		if ( !bed ) {
			cerr << "Error: The requested bed file (" <<bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}

		string bedLine;
		BED bedEntry;                                                                                                                        
		int lineNum = 0;
		
		cout << "<h3>BED Entries from: " << this->bedFile << "</h3>" << endl;

		while (getline(bed, bedLine)) {

			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {
				vector<string> bedFields;
				Tokenize(bedLine,bedFields);

				lineNum++;
			
				if (this->bed->parseBedLine(bedEntry, bedFields, lineNum)) {
					bedEntry.count = 0; 
					WriteURL(bedEntry, base);
				}
			}
		}
		cout << "</table>" << endl;
		cout << "</p>" << endl;
		cout << "\t</body>" << endl <<"</html>" << endl;  
	}
	// Case 2: STDIN.
	else {
		string bedLine;
		BED bedEntry;                                                                                                                        
		int lineNum = 0;

		cout << "<h3>BED Entries from: stdin </h3>" << endl;

		while (getline(cin, bedLine)) {

			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {
				vector<string> bedFields;
				Tokenize(bedLine,bedFields);

				lineNum++;
				if (this->bed->parseBedLine(bedEntry, bedFields, lineNum)) {
					bedEntry.count = 0;
					WriteURL(bedEntry, base);
				}
			}
		}
		cout << "</table>" << endl;
		cout << "</p>" << endl;
		cout << "\t</body>" << endl <<"</html>" << endl;
	}
}
