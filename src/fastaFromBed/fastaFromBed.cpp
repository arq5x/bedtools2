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

Bed2Fa::Bed2Fa(bool &useName, string &dbFile, string &bedFile, string &fastaOutFile, 
	bool &useFasta, bool &useStrand) {

	if (useName) {
		_useName = true;
	}
	
	_dbFile = dbFile;
	_bedFile = bedFile;
	_fastaOutFile = fastaOutFile;
	_useFasta = useFasta;
	_useStrand = useStrand;
		
	_bed = new BedFile(_bedFile);
	
	ExtractDNA();
}


Bed2Fa::~Bed2Fa(void) {
}


//******************************************************************************
// ExtractDNA
//******************************************************************************
void Bed2Fa::ExtractDNA() {

	/* Make sure that we can oen all of the files successfully*/
	
	// open the fasta database for reading
	ifstream faDb(_dbFile.c_str(), ios::in);
	if ( !faDb ) {
		cerr << "Error: The requested fasta database file (" << _dbFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	
	// open the fasta database for reading
	ofstream faOut(_fastaOutFile.c_str(), ios::out);
	if ( !faOut ) {
		cerr << "Error: The requested fasta output file (" << _fastaOutFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}	
	
	// load the BED file into an unbinned map.
	_bed->loadBedFileIntoMapNoBin();

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

				vector<BED> bedList = _bed->bedMapNoBin[currChrom];

				// loop through each BED entry for this chrom and 
				// print the sequence
				for (unsigned int i = 0; i < bedList.size(); i++) {
					
					unsigned int start = bedList[i].start;
					unsigned int end = bedList[i].end;
					
					if ( (start <= currDNA.size()) && (end <= currDNA.size()) ) {
							
						string dna = currDNA.substr(bedList[i].start, ((bedList[i].end - bedList[i].start)));
						// revcomp if necessary.  Thanks to Thomas Doktor.
						if ((_useStrand == true) && (bedList[i].strand == "-")) {
							reverseComplement(dna);
						}
						
						if (!(_useName)) {
							if (_useFasta == true) {
								if (_useStrand == true) {
							    	faOut << ">" << currChrom << ":" << bedList[i].start << "-" 
								          << bedList[i].end   << "(" << bedList[i].strand << ")" << endl << dna << endl;
								}
								else {
							    	faOut << ">" << currChrom << ":" << bedList[i].start << "-" 
								          << bedList[i].end << endl << dna << endl;	
								}
							}
							else {
								if (_useStrand == true) {
									faOut << currChrom << ":" << bedList[i].start << "-" 
									      << bedList[i].end << "(" << bedList[i].strand << ")"
										  << "\t" << dna << endl;								
								}
								else {
									faOut << currChrom << ":" << bedList[i].start << "-" << bedList[i].end 
										  << "\t" << dna << endl;
								}
							}
					  	}
					  	else {
							if (_useFasta == true) {
					    		faOut << ">" << bedList[i].name << endl << dna << endl;
					  		}
							else {
								faOut << bedList[i].name << "\t" << dna << endl;
							}
						}
					}
					else cerr << "Feature (" << bedList[i].chrom << ":" << start << "-" << end << ") beyond " 
						<< currChrom << " size (" << currDNA.size() << " bp).  Skipping." << endl;
				}
				currDNA = "";	
			}
			currChrom = fastaDbLine.substr(1, fastaDbLine.find_first_of(" ")-1);
		}
	}
	
	// process the last chromosome in the fasta file.
	if (currDNA.size() > 0) {
		
		vector<BED> bedList = _bed->bedMapNoBin[currChrom];

		// loop through each BED entry for this chrom and 
		// print the sequence.
		for (unsigned int i = 0; i < bedList.size(); i++) {
			
			unsigned int start = bedList[i].start;
			unsigned int end = bedList[i].end;

			if ( (start <= currDNA.size()) && (end <= currDNA.size()) ) {			

				string dna = currDNA.substr(bedList[i].start, ((bedList[i].end - bedList[i].start)));
				// revcomp if necessary.  Thanks to Thomas Doktor.
				if ((_useStrand == true) && (bedList[i].strand == "-")) {
					reverseComplement(dna);
				}
				
				if (!(_useName)) {
					if (_useFasta == true) {
						if (_useStrand == true) {
					    	faOut << ">" << currChrom << ":" << bedList[i].start << "-" 
						          << bedList[i].end   << "(" << bedList[i].strand << ")" << endl << dna << endl;
						}
						else {
					    	faOut << ">" << currChrom << ":" << bedList[i].start << "-" 
						          << bedList[i].end << endl << dna << endl;	
						}
					}
					else {
						if (_useStrand == true) {
							faOut << currChrom << ":" << bedList[i].start << "-" 
							      << bedList[i].end << "(" << bedList[i].strand << ")"
								  << "\t" << dna << endl;								
						}
						else {
							faOut << currChrom << ":" << bedList[i].start << "-" << bedList[i].end 
								  << "\t" << dna << endl;
						}
					}
			  	}
			  	else {
					if (_useFasta == true) {
			    		faOut << ">" << bedList[i].name << endl << dna << endl;
			  		}
					else {
						faOut << bedList[i].name << "\t" << dna << endl;
					}
				}
			}
			else cerr << "Feature (" << bedList[i].chrom << ":" << start << "-" << end << ") beyond " 
				<< currChrom << " size (" << currDNA.size() << " bp).  Skipping." << endl;

		}
		currDNA = "";	
	}
}



