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


Bed2Fa::Bed2Fa(bool useName, const string &dbFile, const string &bedFile,
    const string &fastaOutFile, bool useFasta, bool useStrand) {

    _useName      = useName;
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
void Bed2Fa::ReportDNA(const BED &bed, string &dna) {

    // revcomp if necessary.  Thanks to Thomas Doktor.
    if ((_useStrand == true) && (bed.strand == "-"))
        reverseComplement(dna);

    if (!(_useName)) {
        if (_useFasta == true) {
            if (_useStrand == true)
                *_faOut << ">" << bed.chrom << ":" << bed.start << "-" << bed.end   << "(" << bed.strand << ")" << endl << dna << endl;
            else
                *_faOut << ">" << bed.chrom << ":" << bed.start << "-" << bed.end << endl << dna << endl;
        }
        else {
            if (_useStrand == true)
                *_faOut << bed.chrom << ":" << bed.start << "-" << bed.end << "(" << bed.strand << ")" << "\t" << dna << endl;
            else
                *_faOut << bed.chrom << ":" << bed.start << "-" << bed.end << "\t" << dna << endl;
        }
    }
    else {
        if (_useFasta == true)
            *_faOut << ">" << bed.name << endl << dna << endl;
        else
            *_faOut << bed.name << "\t" << dna << endl;
    }
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

    // open and memory-map genome file
    FastaReference *fr = new FastaReference;
    bool memmap = true;
    fr->open(_dbFile, memmap);

    BED bed, nullBed;
    string sequence;

    _bed->Open();
    while (_bed->GetNextBed(bed)) {
        if (_bed->_status == BED_VALID) {
            // make sure we are extracting >= 1 bp
            if (bed.zeroLength == false) {
    
                size_t seqLength = fr->sequenceLength(bed.chrom);
                // seqLength > 0 means chrom was found in index.
                // seqLength == 0 otherwise.
                if (seqLength) {
                    // make sure this feature will not exceed the end of the chromosome.
                    if ( (bed.start <= seqLength) && (bed.end <= seqLength) ) 
                    {
                        int length = bed.end - bed.start;
                        sequence = fr->getSubSequence(bed.chrom, bed.start, length);
                        ReportDNA(bed, sequence);
                    }
                    else
                    {
                        cerr << "Feature (" << bed.chrom << ":" << bed.start << "-" << bed.end << ") beyond the length of "
                            << bed.chrom << " size (" << seqLength << " bp).  Skipping." << endl;
                    }
                }
                else
                {
                    cerr << "WARNING. chromosome (" << bed.chrom << 
                            ") was not found in the FASTA file. Skipping."<< endl;
                }
            }
            // handle zeroLength 
            else {
                cerr << "Feature (" << bed.chrom << ":" << bed.start+1 << "-" << bed.end-1 << ") has length = 0, Skipping." << endl;
            }
            bed = nullBed;
        }
    }
    _bed->Close();
}



