/*
***************************************************************************
bamToFastq.h (c) 2009 Aaron Quinlan

Hall Lab
Department of Biochemistry and Molecular Genetics
University of Virginia

All rights reserved.

Filters BAM alignments based upon user-defined criteria.
***************************************************************************
*/
#include "api/BamAux.h"
#include "api/BamReader.h"
using namespace BamTools;

#include "sequenceUtils.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <map>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BamToFastq {

public:

    // constructor 
    BamToFastq(string bamFile, string fastq1, string fastq2, bool useMateTags, bool pairedEnd);

    // destructor
    ~BamToFastq(void);

        
private:
        
    void SingleFastq();
    void PairedFastq();
    void PairedFastqUseTags();

    string _bamFile;

    BamAlignment _end1;
    BamAlignment _end2;

    string _fastq1, _fastq2;    // the names of the fastq output files
    bool _useMateTags;          // whether or not the mate sequence should be 
                                // extracted from the R2 BAM tag.
    bool _pairedEnd;
};
