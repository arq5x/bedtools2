/*
  ***************************************************************************
   bamToFastqMain.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.
 ***************************************************************************
*/

#define PROGRAM_NAME "bamToFastq"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "bamToFastq.h"
#include "version.h"
using namespace std;


// function declarations
void bamtofastq_help(void);
    

int bamtofastq_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    bool haveInBam     = false;
    bool haveFastq1    = false;
    bool haveFastq2    = false;
    bool useMateTags   = false;
    bool pairedEnd     = false;
        
    // input files
    string inBamFile;
    
    //output files
    string fastq1, fastq2;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bamtofastq_help();


    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if (PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveInBam = true;
                inBamFile = argv[i + 1];
                i++;
            }
        }
        else if (PARAMETER_CHECK("-fq", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveFastq1 = true;
                fastq1 = argv[i + 1];
                i++;
            }
        }
        else if (PARAMETER_CHECK("-fq2", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveFastq2 = true;
                pairedEnd  = true;
                fastq2 = argv[i + 1];
                i++;
            }
        }
        else if (PARAMETER_CHECK("-tags", 5, parameterLength)) {
            useMateTags = true;
        }                  
        else {
          cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }       
    }

    // make sure we have a BAM file
    if (!haveInBam) {
      cerr << endl << "*****" << endl << "*****ERROR: Need -bam. " << endl << "*****" << endl;
      showHelp = true;
    }
    // make sure we have an end1 FASTQ file
    if (!haveFastq1) {
      cerr << endl << "*****" << endl << "*****ERROR: Need -fq. " << endl << "*****" << endl;
      showHelp = true;
    }
    
    // let 'er rip.
    if (!showHelp) {
        BamToFastq b2fq(inBamFile, fastq1, fastq2, useMateTags, pairedEnd);
    }
    else {
        bamtofastq_help();
    }
    return 0;
}


void bamtofastq_help(void) {

    cerr << "\nTool:    bedtools bamtofastq (aka bamToFastq)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Convert BAM alignments to FASTQ files. " << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <BAM> -fq <FQ> " << endl << endl;

    cerr << "Options:" << endl;
    cerr << "\t-fq2\tFASTQ for second end.  Used if BAM contains paired-end data." << endl;
    cerr << "\t\tBAM should be sorted by query name is creating paired FASTQ." << endl << endl;
    
    cerr << "\t-tags\tCreate FASTQ based on the mate info" << endl;
    cerr << "\t\tin the BAM R2 and Q2 tags." << endl << endl;
    
    cerr << "Tips: " << endl;
    cerr << "\tIf you want to create a single, interleaved FASTQ file " << endl;
    cerr << "\tfor paired-end data, you can just write both to /dev/stdout:" << endl << endl;
    cerr << "\tbedtools bamtofastq -i x.bam -fq /dev/stdout -fq2 /dev/stdout > x.ilv.fq" << endl << endl;
    exit(1);    
}

