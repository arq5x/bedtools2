/*****************************************************************************
  linksBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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
		else if ((this->bed->bedType == 6) || (this->bed->bedType == 9) || (this->bed->bedType == 12)) {
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


void BedLinks::LinksBed(istream &bedInput) {


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
	
	string bedLine;
	BED bedEntry;                                                                                                                        
	int lineNum = 0;
	vector<string> bedFields;
	bedFields.reserve(12);
	
	cout << "<h3>BED Entries from: stdin </h3>" << endl;

	while (getline(bedInput, bedLine)) {

		Tokenize(bedLine,bedFields);
		lineNum++;
			
		if (this->bed->parseLine(bedEntry, bedFields, lineNum)) {
			bedEntry.count = 0;
			WriteURL(bedEntry, base);
		}
		bedFields.clear();
	}
	cout << "</table>" << endl;
	cout << "</p>" << endl;
	cout << "\t</body>" << endl <<"</html>" << endl;
}


void BedLinks::DetermineBedInput() {

	if (this->bedFile != "stdin") {   // process a file
		ifstream beds(this->bedFile.c_str(), ios::in);
		if ( !beds ) {
			cerr << "Error: The requested bed file (" << this->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		LinksBed(beds);
	}
	else {   // process stdin
		LinksBed(cin);		
	}
}
