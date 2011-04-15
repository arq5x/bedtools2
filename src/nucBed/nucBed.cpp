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


NucBed::NucBed(string &dbFile, string &bedFile, bool printSeq) {

    _dbFile       = dbFile;
    _bedFile      = bedFile;
    _printSeq     = printSeq;

    _bed = new BedFile(_bedFile);

    // Compute the DNA content in each BED/GFF/VCF interval
    ProfileDNA();
}


NucBed::~NucBed(void) 
{}


void NucBed::ReportDnaProfile(const BED& bed, const string &sequence, int seqLength)
{
    int a,c,g,t,n,other;
    a = c = g = t = n = other = 0;
    
    getDnaContent(sequence,a,c,g,t,n,other);
    
    // report the original interval
    _bed->reportBedTab(bed);
    // report AT and GC content
    printf("%f\t%f\t",(float)(a+t)/seqLength, (float)(c+g)/seqLength);
    // report raw nucleotide counts
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d",a,c,g,t,n,other,seqLength);
    // add the original sequence if requested.
    if (_printSeq) {
        printf("\t%s\n",sequence.c_str());
    }
    else {
        printf("\n");
    }
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
    fr.open(_dbFile, memmap);

    BED bed, nullBed;
    int lineNum = 0;
    BedLineStatus bedStatus;
    string sequence;

    _bed->Open();
    while ((bedStatus = _bed->GetNextBed(bed, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            // make sure we are extracting >= 1 bp
            if (bed.zeroLength == false) {
                size_t seqLength = fr.sequenceLength(bed.chrom);
                // make sure this feature will not exceed the end of the chromosome.
                if ( (bed.start <= seqLength) && (bed.end <= seqLength) ) 
                {
                    // grab the dna at this interval
                    int length = bed.end - bed.start;
                    // report the sequence's content
                    ReportDnaProfile(bed, fr.getSubSequence(bed.chrom, bed.start, length), length);
                    bed = nullBed;
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
            bed = nullBed;
        }
    }
    _bed->Close();
}



