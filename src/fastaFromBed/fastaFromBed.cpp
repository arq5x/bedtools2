/*****************************************************************************
  fastaFromBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "fastaFromBed.h"

Bed2Fa::Bed2Fa(bool &useName, string &dbFile, string &bedFile, string &fastaOutFile) {

	if (useName) {
		this->useName = true;
	}
	
	this->dbFile = dbFile;
	this->bedFile = bedFile;
	this->fastaOutFile = fastaOutFile;
	
	this->bed = new BedFile(this->bedFile);
}


Bed2Fa::~Bed2Fa(void) {
}


//******************************************************************************
// ExtractDNA
//******************************************************************************
void Bed2Fa::ExtractDNA() {

	/* Make sure that we can oen all of the files successfully*/
	
	// open the fasta database for reading
	ifstream faDb(this->dbFile.c_str(), ios::in);
	if ( !faDb ) {
		cerr << "Error: The requested fasta database file (" << this->dbFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	
	// open the fasta database for reading
	ofstream faOut(this->fastaOutFile.c_str(), ios::out);
	if ( !faOut ) {
		cerr << "Error: The requested fasta output file (" << this->fastaOutFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}	
	
	// load the BED file into an unbinned map.
	bed->loadBedFileIntoMapNoBin();

	//Read the fastaDb chromosome by chromosome
	string fastaDbLine;
	string currChrom;
	string currDNA = "";
	currDNA.reserve(500000000);
	
	while (getline(faDb,fastaDbLine)) {
		
		if (fastaDbLine.find(">",0) != 0 ) {
			currDNA += fastaDbLine;
		}
		else {
			if (currDNA.size() > 0) {

				vector<BED> bedList = bed->bedMapNoBin[currChrom];

				// loop through each BED entry for this chrom and 
				// print the sequence
				for (unsigned int i = 0; i < bedList.size(); i++) {
						
					string dna = currDNA.substr(bedList[i].start, ((bedList[i].end - bedList[i].start)));
				
					if (!(this->useName)) {
				    	faOut << ">" << currChrom << ":"
					  		<< bedList[i].start << "-" << bedList[i].end << endl << dna << endl;
				  	}
				  	else {
				    	faOut << ">" << bedList[i].name << endl << dna << endl;
				  	}
				}
				currDNA = "";	
			}
			currChrom = fastaDbLine.substr(1, fastaDbLine.find_first_of(" ")-1);
		}
	}
	
	// process the last chromosome in the fasta file.
	if (currDNA.size() > 0) {
		
		vector<BED> bedList = bed->bedMapNoBin[currChrom];

		// loop through each BED entry for this chrom and 
		// print the sequence.
		for (unsigned int i = 0; i < bedList.size(); i++) {
				
			string dna = currDNA.substr(bedList[i].start, ((bedList[i].end - bedList[i].start)));
		
			if (!(this->useName)) {
		    	faOut << ">" << currChrom << ":"
			  		<< bedList[i].start << "-" << bedList[i].end << endl << dna << endl;
		  	}
		  	else {
		    	faOut << ">" << bedList[i].name << endl << dna << endl;
		  	}
		}
		currDNA = "";	
	}
}



