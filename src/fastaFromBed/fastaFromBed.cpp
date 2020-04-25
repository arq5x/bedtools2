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
#include "bedFile.h"


Bed2Fa::Bed2Fa(const string &dbFile, 
               const string &bedFile, const string &fastaOutFile,
               bool useFasta, bool useStrand, 
               bool useBlocks, bool useFullHeader,
               bool useBedOut, bool useName, 
               bool useNamePlus, bool useNameOnly,
               bool isRNA) :
    _dbFile(dbFile),
    _bedFile(bedFile),
    _fastaOutFile(fastaOutFile),
    _useFasta(useFasta),
    _useStrand(useStrand),
    _useBlocks(useBlocks),
    _useFullHeader(useFullHeader),
    _useBedOut(useBedOut),
    _useName(useName),
    _useNamePlus(useNamePlus),
    _useNameOnly(useNameOnly),
    _isRNA(isRNA)
{
    _bed = new BedFile(_bedFile);

      // Figure out what the output file should be.
    if (fastaOutFile == "stdout" || fastaOutFile == "-") {
        _faOut = &cout;
    }
    else {
        // Make sure we can open the file.
        ofstream fa(fastaOutFile.c_str(), ios::out);
        if ( !fa ) {
            cerr << "Error: The requested fasta output file (" 
                 << fastaOutFile << ") could not be opened. Exiting!" 
                 << endl;
            exit (1);
        }
        else {
            fa.close();
            _faOut = new ofstream(fastaOutFile.c_str(), ios::out);
        }
    }
    // Extract the requested intervals from the FASTA input file.
    ExtractSeq();
}


Bed2Fa::~Bed2Fa(void) {
}


//******************************************************************************
// ReportSeq
//******************************************************************************
void Bed2Fa::ReportSeq(const BED &bed, string &seq) {

    // revcomp if necessary.  Thanks to Thomas Doktor.
    if ((_useStrand == true) && (bed.strand == "-"))
        reverseComplement(seq, _isRNA);

    if (_useBedOut) {
        _bed->reportBedTab(bed);
        *_faOut  << seq << endl;
    }
    else {
        ostringstream header;
        if (_useName || _useNamePlus)
        {
            header << bed.name << "::" << bed.chrom << ":" 
                   << bed.start << "-" << bed.end;
        }
        else if (_useNameOnly)
        {
            header << bed.name;
        }
        else 
        {
            header << bed.chrom << ":" 
                   << bed.start << "-" << bed.end;
        }
        
        if (_useStrand)
        {
            header << "(" << bed.strand << ")";
        }
        
        if (_useFasta)
            *_faOut << ">" << header.str() << endl << seq << endl;
        else
            *_faOut << header.str() << "\t" << seq << endl;
    }
}

//******************************************************************************
// Extract sequence
//******************************************************************************
void Bed2Fa::ExtractSeq() {

    /* Make sure that we can open all of the files successfully*/

    // open the fasta database for reading
    ifstream faDb(_dbFile.c_str(), ios::in);
    if ( !faDb ) {
        cerr << "Error: The requested fasta database file (" 
             << _dbFile << ") could not be opened. Exiting!" 
             << endl;
        exit (1);
    }

    // open and memory-map genome file
    FastaReference *fr = new FastaReference;
    bool memmap = true;
    fr->open(_dbFile, memmap, _useFullHeader);

    BED bed, nullBed;
    string sequence;

    _bed->Open();
    while (_bed->GetNextBed(bed)) {
        if (_bed->_status == BED_VALID) {
            // make sure we are extracting >= 1 bp
            if (bed.zeroLength == false) {
    
                CHRPOS seqLength = fr->sequenceLength(bed.chrom);
                // seqLength > 0 means chrom was found in index.
                // seqLength == 0 otherwise.
                if (seqLength) {
                    // make sure this feature will not exceed 
                    // the end of the chromosome.
                    if ( (bed.start <= seqLength) && (bed.end <= seqLength) ) 
                    {
                        CHRPOS length = bed.end - bed.start;
                        if(_useBlocks){
                            // vec to store the discrete BED "blocks"
                            bedVector bedBlocks;  
                            GetBedBlocks(bed, bedBlocks);
                            sequence.clear();
                            for (int i = 0; i < (int) bedBlocks.size(); ++i) {
                                sequence += fr->getSubSequence(bed.chrom,
                                        bedBlocks[i].start,
                                        bedBlocks[i].end - bedBlocks[i].start);
                            }
                        } else {
                            sequence = \
                               fr->getSubSequence(bed.chrom, bed.start, length);
                        }
                        ReportSeq(bed, sequence);
                    }
                    else
                    {
                        cerr << "Feature (" << bed.chrom << ":" 
                             << bed.start << "-" << bed.end 
                            << ") beyond the length of "
                            << bed.chrom 
                            << " size (" << seqLength << " bp).  Skipping." 
                            << endl;
                    }
                }
                else
                {
                    cerr << "WARNING. chromosome (" 
                         << bed.chrom 
                         << ") was not found in the FASTA file. Skipping."
                         << endl;
                }
            }
            // handle zeroLength 
            else {
                cerr << "Feature (" << bed.chrom << ":" 
                     << bed.start+1 << "-" << bed.end-1 
                     << ") has length = 0, Skipping." 
                     << endl;
            }
            bed = nullBed;
        }
    }
    _bed->Close();
}



