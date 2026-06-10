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
               bool useNameKey, const string &nameKey,
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
    _useNameKey(useNameKey),
    _nameKey(nameKey),
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
        if (_useNameKey)
        {
            string attrVal;
            if (getGffAttribute(bed, _nameKey, attrVal)) {
                header << attrVal;
            }
            else {
                cerr << "WARNING: attribute '" << _nameKey
                     << "' not found for feature at " << bed.chrom << ":"
                     << bed.start << "-" << bed.end
                     << "; using feature name." << endl;
                header << bed.name;
            }
        }
        else if (_useName || _useNamePlus)
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


bool Bed2Fa::getGffAttribute(const BED &bed, const string &key, string &value) {
    // Only GFF records carry a column-9 attribute string.
    if (!_bed->isGff() || bed.fields.size() < 9)
        return false;

    const string &attrs = bed.fields[8];
    size_t pos = 0;
    while (pos < attrs.size()) {
        size_t semi = attrs.find(';', pos);
        string token = attrs.substr(pos, (semi == string::npos)
                                         ? string::npos : semi - pos);
        // trim leading/trailing whitespace
        size_t b = token.find_first_not_of(" \t");
        size_t e = token.find_last_not_of(" \t");
        if (b != string::npos) {
            token = token.substr(b, e - b + 1);
            size_t eq = token.find('=');
            if (eq != string::npos && token.substr(0, eq) == key) {
                string v = token.substr(eq + 1);
                // strip one pair of surrounding double-quotes
                if (v.size() >= 2 && v.front() == '"' && v.back() == '"')
                    v = v.substr(1, v.size() - 2);
                value = v;
                return true;
            }
        }
        if (semi == string::npos) break;
        pos = semi + 1;
    }
    return false;
}



