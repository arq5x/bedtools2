/*****************************************************************************
  nucBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "nucBed.h"


NucBed::NucBed(string &dbFile, string &bedFile, bool printSeq, 
               bool hasPattern, const string &pattern, 
               bool forceStrand, bool ignoreCase, bool useFullHeader) 
{

    _dbFile        = dbFile;
    _bedFile       = bedFile;
    _printSeq      = printSeq;
    _hasPattern    = hasPattern;
    _pattern       = pattern;
    _forceStrand   = forceStrand;
    _ignoreCase    = ignoreCase;
    _useFullHeader = useFullHeader;
    _bed = new BedFile(_bedFile);

    // Compute the DNA content in each BED/GFF/VCF interval
    ProfileDNA();
}


NucBed::~NucBed(void) 
{}


void NucBed::ReportDnaProfile(const BED& bed, const string &sequence, CHRPOS seqLength)
{
    CHRPOS a,c,g,t,n,other,userPatternCount;
    a = c = g = t = n = other = userPatternCount = 0;
    
    getDnaContent(sequence,a,c,g,t,n,other);
    
    if (_hasPattern)
        userPatternCount = countPattern(sequence, _pattern, _ignoreCase);
    
    
    // report the original interval
    _bed->reportBedTab(bed);
    // report AT and GC content
    printf("%f\t%f\t",(float)(a+t)/seqLength, (float)(c+g)/seqLength);
    // report raw nucleotide counts
    printf("%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%" PRId_CHRPOS,a,c,g,t,n,other,seqLength);
    // add the original sequence if requested.

    if (_printSeq)
        printf("\t%s",sequence.c_str());
    if (_hasPattern)
        printf("\t%" PRId_CHRPOS,userPatternCount);
    printf("\n");

}


void NucBed::PrintHeader(void) {
    printf("#");
    
    int numOrigColumns = (int) _bed->bedType;
    for (int i = 1; i <= numOrigColumns; ++i) {
        printf("%d_usercol\t", i);
    }
    printf("%d_pct_at\t", numOrigColumns + 1);
    printf("%d_pct_gc\t", numOrigColumns + 2);
    printf("%d_num_A\t", numOrigColumns + 3);
    printf("%d_num_C\t", numOrigColumns + 4);
    printf("%d_num_G\t", numOrigColumns + 5);
    printf("%d_num_T\t", numOrigColumns + 6);
    printf("%d_num_N\t", numOrigColumns + 7);
    printf("%d_num_oth\t", numOrigColumns + 8);
    printf("%d_seq_len", numOrigColumns + 9);
    
    if (_printSeq)
        printf("\t%d_seq", numOrigColumns + 10);
    if (_hasPattern && !_printSeq)
        printf("\t%d_user_patt_count", numOrigColumns + 10);
    else if (_hasPattern && _printSeq)
        printf("\t%d_user_patt_count", numOrigColumns + 11);
    printf("\n");

}


//******************************************************************************
// ExtractDNA
//******************************************************************************
void NucBed::ProfileDNA() {

    /* Make sure that we can oen all of the files successfully*/

    // open the fasta database for reading
    ifstream faDb(_dbFile.c_str(), ios::in);
    if ( !faDb ) {
        cerr << "Error: The requested fasta database file (" << _dbFile << ") could not be opened. Exiting!" << endl;
        exit (1);
    }

    // open and memory-map genome file
    FastaReference fr;
    bool memmap = true;    
    fr.open(_dbFile, memmap, _useFullHeader);

    bool headerReported = false;
    BED bed;
    string sequence;

    _bed->Open();
    while (_bed->GetNextBed(bed)) {
        if (_bed->_status == BED_VALID) {
            if (headerReported == false) {
                PrintHeader();
                headerReported = true;
            }
            // make sure we are extracting >= 1 bp
            if (bed.zeroLength == false) {
                CHRPOS seqLength = fr.sequenceLength(bed.chrom);
                // make sure this feature will not exceed the end of the chromosome.
                if ( (bed.start <= seqLength) && (bed.end <= seqLength) ) 
                {
                    // grab the dna at this interval
                    CHRPOS length = bed.end - bed.start;
                    // report the sequence's content
                    string dna = fr.getSubSequence(bed.chrom, bed.start, length);
                    // rev comp si necessaire
                    if ((_forceStrand == true) && (bed.strand == "-"))
                        reverseComplement(dna);
                    ReportDnaProfile(bed, dna, length);
                }
                else
                {
                    cerr << "Feature (" << bed.chrom << ":" << bed.start << "-" << bed.end << ") beyond the length of "
                        << bed.chrom << " size (" << seqLength << " bp).  Skipping." << endl;
                }
            }
            // handle zeroLength 
            else {
                cerr << "Feature (" << bed.chrom << ":" << bed.start+1 << "-" << bed.end-1 << ") has length = 0, Skipping." << endl;
            }
        }
    }
    _bed->Close();
}



