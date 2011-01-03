/*****************************************************************************
  fastaFromBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "fastaFromBed.h"


Bed2Fa::Bed2Fa(bool &useName, string &dbFile, string &bedFile,
    string &fastaOutFile, bool &useFasta, bool &useStrand) {

    if (useName) {
        _useName = true;
    }

    _dbFile       = dbFile;
    _bedFile      = bedFile;
    _fastaOutFile = fastaOutFile;
    _useFasta     = useFasta;
    _useStrand    = useStrand;

    _bed = new BedFile(_bedFile);

    // Figure out what the output file should be.
    if (fastaOutFile == "stdout") {
        _faOut = &cout;
    }
    else {
        // Make sure we can open the file.
        ofstream fa(fastaOutFile.c_str(), ios::out);
        if ( !fa ) {
            cerr << "Error: The requested fasta output file (" << fastaOutFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            fa.close();
            _faOut = new ofstream(fastaOutFile.c_str(), ios::out);
        }
    }

    // Extract the requested intervals from the FASTA input file.
    ExtractDNA();
}


Bed2Fa::~Bed2Fa(void) {
}


//******************************************************************************
// ReportDNA
//******************************************************************************
void Bed2Fa::ReportDNA(const BED &bed, const string &currDNA, const string &currChrom) {

    if ( (bed.start <= currDNA.size()) && (bed.end <= currDNA.size()) ) {

        string dna = currDNA.substr(bed.start, ((bed.end - bed.start)));
        // revcomp if necessary.  Thanks to Thomas Doktor.
        if ((_useStrand == true) && (bed.strand == "-"))
            reverseComplement(dna);

        if (!(_useName)) {
            if (_useFasta == true) {
                if (_useStrand == true)
                    *_faOut << ">" << currChrom << ":" << bed.start << "-" << bed.end   << "(" << bed.strand << ")" << endl << dna << endl;
                else
                    *_faOut << ">" << currChrom << ":" << bed.start << "-" << bed.end << endl << dna << endl;
            }
            else {
                if (_useStrand == true)
                    *_faOut << currChrom << ":" << bed.start << "-" << bed.end << "(" << bed.strand << ")" << "\t" << dna << endl;
                else
                    *_faOut << currChrom << ":" << bed.start << "-" << bed.end << "\t" << dna << endl;
            }
        }
        else {
            if (_useFasta == true)
                *_faOut << ">" << bed.name << endl << dna << endl;
            else
                *_faOut << bed.name << "\t" << dna << endl;
        }
    }
    else cerr << "Feature (" << bed.chrom << ":" << bed.start << "-" << bed.end << ") beyond "
        << currChrom << " size (" << currDNA.size() << " bp).  Skipping." << endl;
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

                vector<BED>::const_iterator bedItr = _bed->bedMapNoBin[currChrom].begin();
                vector<BED>::const_iterator bedEnd = _bed->bedMapNoBin[currChrom].end();
                // loop through each BED entry for this chrom and print the sequence
                for (; bedItr != bedEnd; ++bedItr) {
                    ReportDNA(*bedItr, currDNA, currChrom);
                }
                currDNA = "";
            }
            currChrom = fastaDbLine.substr(1, fastaDbLine.find_first_of(" ")-1);
        }
    }

    // process the last chromosome in the fasta file.
    if (currDNA.size() > 0) {
        vector<BED>::const_iterator bedItr = _bed->bedMapNoBin[currChrom].begin();
        vector<BED>::const_iterator bedEnd = _bed->bedMapNoBin[currChrom].end();
        // loop through each BED entry for this chrom and print the sequence
        for (; bedItr != bedEnd; ++bedItr) {
            ReportDNA(*bedItr, currDNA, currChrom);
        }
        currDNA = "";
    }
}



