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
    _fq = new ofstream(_fastq1.c_str(), ios::out);
    if ( !*_fq ) {
        cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the BAM file
    BamReader reader;
    reader.Open(_bamFile);
    if (!reader.Open(_bamFile)) {
        cerr << "Failed to open BAM file " << _bamFile << endl;
        exit(1);
    }
    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
        // extract the sequence and qualities for the BAM "query"
        string seq  = bam.QueryBases;
        string qual = bam.Qualities;
        if (bam.IsReverseStrand() == true) {
            reverseComplement(seq);
            reverseSequence(qual);
        }

        *_fq << "@" << bam.Name << endl;
        *_fq << seq << endl;
        *_fq << "+" << endl;
        *_fq << qual << endl;
    }
}

void BamToFastq::PairedFastq() {
    // open the 1st fastq file for writing
    _fq1 = new ofstream(_fastq1.c_str(), ios::out);
    if ( !*_fq1 ) {
        cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the 2nd fastq file for writing
    _fq2 = new ofstream(_fastq2.c_str(), ios::out);
    if ( !*_fq2 ) {
        cerr << "Error: The second fastq file (" << _fastq2 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the BAM file
    BamReader reader;
    if (!reader.Open(_bamFile)) {
        cerr << "Failed to open BAM file " << _bamFile << endl;
        exit(1);
    }
    // rip through the BAM file and convert each mapped entry to BEDPE
    vector<BamAlignment> alignments;
    string prevName = "";
    string currName = "";
    BamAlignment curr;
    while (reader.GetNextAlignment(curr)) {
        currName = curr.Name;
        if ((currName != prevName) && (prevName != "")) 
        {
            WritePairs(alignments);
            alignments.clear();
            alignments.push_back(curr);
        }
        else {
            if (curr.IsPaired())
            {
                alignments.push_back(curr);
            }
        }
        prevName = currName;
    }
    if (alignments.size() != 0) WritePairs(alignments);
    reader.Close();
}

void BamToFastq::WritePairs(vector<BamAlignment> &alignments) 
{
    string seq1, seq2;
    string qual1, qual2;
    string queryName;
    if (alignments.size() == 1)
    {
        cerr << "*****WARNING: Query " << alignments[0].Name
             << " is marked as paired, but its mate does not occur"
             << " next to it in your BAM file.  Skipping. " << endl;
        return;
    }
    else if (alignments.size() == 2)
    {
        queryName = alignments[0].Name;
        BamAlignment curr = alignments[0];
        if (alignments[0].IsFirstMate() && alignments[1].IsSecondMate())
        {
            seq1  = alignments[0].QueryBases;
            qual1 = alignments[0].Qualities;
            seq2  = alignments[1].QueryBases;
            qual2 = alignments[1].Qualities;
            if (alignments[0].IsReverseStrand() == true) {
                reverseComplement(seq1);
                reverseSequence(qual1);
            }
            if (alignments[1].IsReverseStrand() == true) {
                reverseComplement(seq2);
                reverseSequence(qual2);
            }
        }
        else if (alignments[0].IsSecondMate() && alignments[1].IsFirstMate())
        {
            seq1  = alignments[1].QueryBases;
            qual1 = alignments[1].Qualities;
            seq2  = alignments[0].QueryBases;
            qual2 = alignments[0].Qualities;
            if (alignments[1].IsReverseStrand() == true) {
                reverseComplement(seq1);
                reverseSequence(qual1);
            }
            if (alignments[0].IsReverseStrand() == true) {
                reverseComplement(seq2);
                reverseSequence(qual2);
            }
        }
        else {
            cerr << "*****WARNING: Query " << alignments[0].Name
                 << " is marked as paired, but first and/or second mates"
                 << " were not found in your BAM file.  Skipping. " << endl;
            return;
        }
    }
    else if (alignments.size() > 2)
    {
        queryName = alignments[0].Name;
        size_t mateIndex = 0;
        if (alignments[0].IsFirstMate())
        {
            for (size_t i = 1; i < alignments.size(); ++i)
            {
                if (alignments[i].IsSecondMate())
                {
                    mateIndex = i;
                    break;
                }
            }
            if (mateIndex > 0)
            {
                seq1  = alignments[0].QueryBases;
                qual1 = alignments[0].Qualities;
                seq2  = alignments[mateIndex].QueryBases;
                qual2 = alignments[mateIndex].Qualities;
                if (alignments[0].IsReverseStrand() == true) {
                    reverseComplement(seq1);
                    reverseSequence(qual1);
                }
                if (alignments[mateIndex].IsReverseStrand() == true) {
                    reverseComplement(seq2);
                    reverseSequence(qual2);
                }
            }
            else {
                cerr << "*****WARNING: Query " << alignments[0].Name
                    << " is marked as paired, but only was found"
                    << " in your BAM file.  Skipping. " << endl;
                return;
            }
        }
        else if (alignments[0].IsSecondMate())
        {
            for (size_t i = 1; i < alignments.size(); ++i)
            {
                if (alignments[i].IsFirstMate())
                {
                    mateIndex = i;
                    break;
                }
            }
            if (mateIndex > 0)
            {
                seq2  = alignments[0].QueryBases;
                qual2 = alignments[0].Qualities;
                seq1  = alignments[mateIndex].QueryBases;
                qual1 = alignments[mateIndex].Qualities;
                if (alignments[0].IsReverseStrand() == true) {
                    reverseComplement(seq1);
                    reverseSequence(qual1);
                }
                if (alignments[mateIndex].IsReverseStrand() == true) {
                    reverseComplement(seq2);
                    reverseSequence(qual2);
                }
            }
            else {
                cerr << "*****WARNING: Query " << alignments[0].Name
                    << " is marked as paired, but only was found"
                    << " in your BAM file.  Skipping. " << endl;
                return;
            }
        }
        else {
            cerr << "*****WARNING: Query " << alignments[0].Name
                 << " is marked as paired, but first and/or second mates"
                 << " were not found in your BAM file.  Skipping. " << endl;
            return;
        }
    }
    *_fq1 << "@" << queryName << "/1" << endl;
    *_fq1 << seq1 << endl;
    *_fq1 << "+" << endl;
    *_fq1 << qual1 << endl;
    
    *_fq2 << "@" << queryName << "/2" << endl;
    *_fq2 << seq2 << endl;
    *_fq2 << "+" << endl;
    *_fq2 << qual2 << endl;
}

void BamToFastq::PairedFastqUseTags() {

    // open the 1st fastq file for writing
    _fq1 = new ofstream(_fastq1.c_str(), ios::out);
    if ( !*_fq1 ) {
        cerr << "Error: The first fastq file (" << _fastq1 << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    // open the 2nd fastq file for writing
    _fq2 = new ofstream(_fastq2.c_str(), ios::out);
    if ( !*_fq2 ) {
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
                *_fq1 << "@" << bam1.Name << "/1" << endl;
                *_fq1 << seq1 << endl;
                *_fq1 << "+" << endl;
                *_fq1 << qual1 << endl;
                // end2
                *_fq2 << "@" << bam1.Name << "/2" <<endl;
                *_fq2 << mateSequence << endl;
                *_fq2 << "+" << endl;
                *_fq2 << mateQualities << endl;
            }
            else {
                // end 2
                *_fq2 << "@" << bam1.Name << "/2" <<endl;
                *_fq2 << seq1 << endl;
                *_fq2 << "+" << endl;
                *_fq2 << qual1 << endl;
                // end 1
                *_fq1 << "@" << bam1.Name << "/1" <<endl;
                *_fq1 << mateSequence << endl;
                *_fq1 << "+" << endl;
                *_fq1 << mateQualities << endl;
            }
        }
    }
    reader.Close();
}


