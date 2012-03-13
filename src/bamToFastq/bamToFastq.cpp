/*
  ***************************************************************************
   bamToFastq.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.

   Filters BAM alignments based upon user-defined criteria.
 ***************************************************************************
*/

#include "bamToFastq.h"

// constructor
BamToFastq::BamToFastq(string bamFile, string fastq1, string fastq2, bool useMateTags, bool pairedEnd)
: _bamFile(bamFile)
, _fastq1(fastq1)
, _fastq2(fastq2)
, _useMateTags(useMateTags)
, _pairedEnd(pairedEnd)
{
    if (!_pairedEnd)
        SingleFastq();
    else {
        if (!_useMateTags) PairedFastq();
        else PairedFastqUseTags();
    }
}


// destructor
BamToFastq::~BamToFastq(void) {}


void BamToFastq::SingleFastq() {
    // open the 1st fastq file for writing
    ofstream fq(_fastq1.c_str(), ios::out);
    if ( !fq ) {
        cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the BAM file
    BamReader reader;
    reader.Open(_bamFile);
    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
        // extract the sequence and qualities for the BAM "query"
        string seq  = bam.QueryBases;
        string qual = bam.Qualities;
        if (bam.IsReverseStrand() == true) {
            reverseComplement(seq);
            reverseSequence(qual);
        }
        fq << "@" << bam.Name << endl;
        fq << seq << endl;
        fq << "+" << endl;
        fq << qual << endl;
    }
}

void BamToFastq::PairedFastq() {
    // open the 1st fastq file for writing
    ofstream fq1(_fastq1.c_str(), ios::out);
    if ( !fq1 ) {
        cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the 2nd fastq file for writing
    ofstream fq2(_fastq2.c_str(), ios::out);
    if ( !fq2 ) {
        cerr << "Error: The second fastq file (" << _fastq2 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the BAM file
    BamReader reader;
    reader.Open(_bamFile);
    // rip through the BAM file and convert each mapped entry to BEDPE
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        
        reader.GetNextAlignment(bam2);        
        if (bam1.Name != bam2.Name) {
            while (bam1.Name != bam2.Name)
            {
                if (bam1.IsPaired()) 
                {
                    cerr << "*****WARNING: Query " << bam1.Name
                         << " is marked as paired, but it's mate does not occur"
                         << " next to it in your BAM file.  Skipping. " << endl;
                }
                bam1 = bam2;
                reader.GetNextAlignment(bam2);
            }
        }
        else if (bam1.IsPaired() && bam2.IsPaired()) {
            // extract the sequence and qualities for the BAM "query"
            string seq1  = bam1.QueryBases;
            string qual1 = bam1.Qualities;
            string seq2  = bam2.QueryBases;
            string qual2 = bam2.Qualities;
            if (bam1.IsReverseStrand() == true) {
                reverseComplement(seq1);
                reverseSequence(qual1);
            }
            if (bam2.IsReverseStrand() == true) {
                reverseComplement(seq2);
                reverseSequence(qual2);
            }
            fq1 << "@" << bam1.Name << "/1" << endl;
            fq1 << seq1 << endl;
            fq1 << "+" << endl;
            fq1 << qual1 << endl;
            
            fq2 << "@" << bam2.Name << "/2" << endl;
            fq2 << seq2 << endl;
            fq2 << "+" << endl;
            fq2 << qual2 << endl;
        }
    }
    reader.Close();
}


void BamToFastq::PairedFastqUseTags() {

    // open the 1st fastq file for writing
    ofstream fq1(_fastq1.c_str(), ios::out);
    if ( !fq1 ) {
        cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the 2nd fastq file for writing
    ofstream fq2(_fastq2.c_str(), ios::out);
    if ( !fq2 ) {
        cerr << "Error: The second fastq file (" << _fastq2 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }

    // open the BAM file
    BamReader reader;
    reader.Open(_bamFile);
    // rip through the BAM file and convert each mapped entry to BEDPE
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        
        reader.GetNextAlignment(bam2);        
        if (bam1.Name != bam2.Name) {
            while (bam1.Name != bam2.Name)
            {
                if (bam1.IsPaired()) 
                {
                    cerr << "*****WARNING: Query " << bam1.Name
                         << " is marked as paired, but it's mate does not occur"
                         << " next to it in your BAM file.  Skipping. " << endl;
                }
                bam1 = bam2;
                reader.GetNextAlignment(bam2);
            }
        }
        else if (bam1.IsPaired() && bam2.IsPaired()) {
            // assume the R2 and Q2 tags are on the + strand.
            string mateSequence, mateQualities;
            bam1.GetTag("R2", mateSequence);
            bam1.GetTag("Q2", mateQualities);

            string seq1  = bam1.QueryBases;
            string qual1 = bam1.Qualities;
            if (bam1.IsReverseStrand() == true) {
                reverseComplement(seq1);
                reverseSequence(qual1);
            }
            
            // since the info for both ends are contained in each BAM record,
            // we only need to process one of the two records (bam1) in order
            // to produce FASTQ entries for both ends.
            // NOTE: Assumes that R2 and Q2 have already been rev 
            //      and revcomped if necessary
            if (bam1.IsFirstMate() == true) {
                // end1
                fq1 << "@" << bam1.Name << "/1" << endl;
                fq1 << seq1 << endl;
                fq1 << "+" << endl;
                fq1 << qual1 << endl;
                // end2
                fq2 << "@" << bam1.Name << "/2" <<endl;
                fq2 << mateSequence << endl;
                fq2 << "+" << endl;
                fq2 << mateQualities << endl;
            }
            else {
                // end 2
                fq2 << "@" << bam1.Name << "/2" <<endl;
                fq2 << seq1 << endl;
                fq2 << "+" << endl;
                fq2 << qual1 << endl;
                // end 1
                fq1 << "@" << bam1.Name << "/1" <<endl;
                fq1 << mateSequence << endl;
                fq1 << "+" << endl;
                fq1 << mateQualities << endl;
            }
        }
    }
    reader.Close();
}


