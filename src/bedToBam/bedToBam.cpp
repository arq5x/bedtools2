/*****************************************************************************
  bedToBam.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "bedFile.h"
#include "genomeFile.h"
#include "version.h"


#include "api/BamReader.h"
#include "api/BamAux.h"
#include "api/BamWriter.h"
using namespace BamTools;

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools bedtobam"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void bedtobam_help(void);
void ProcessBed(BedFile *bed, GenomeFile *genome, bool isBED12, int mapQual, bool uncompressedBam);
void ConvertBedToBam(const BED &bed, BamAlignment &bam, map<string, int> &chromToId, bool isBED12, int mapQual, int lineNum);
void MakeBamHeader(const string &genomeFile, RefVector &refs, string &header, map<string, int> &chromToInt);
int  bedtobam_reg2bin(int beg, int end);



int bedtobam_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile = "stdin";
    string genomeFile;

    unsigned int mapQual = 255;

    bool haveBed         = true;
    bool haveGenome      = false;
    bool haveMapQual     = false;
    bool isBED12         = false;
    bool uncompressedBam = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bedtobam_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveGenome = true;
                genomeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-mapq", 5, parameterLength)) {
            haveMapQual = true;
            if ((i+1) < argc) {
                mapQual = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bed12", 6, parameterLength)) {
            isBED12 = true;
        }
        else if(PARAMETER_CHECK("-ubam", 5, parameterLength)) {
            uncompressedBam = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (!haveBed ) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i (BED) file. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (!haveGenome ) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -g (genome) file. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (mapQual < 0 || mapQual > 255) {
        cerr << endl << "*****" << endl << "*****ERROR: MAPQ must be in range [0,255]. " << endl << "*****" << endl;
        showHelp = true;
    }


    if (!showHelp) {
        BedFile *bed       = new BedFile(bedFile);
        GenomeFile *genome = new GenomeFile(genomeFile);

        ProcessBed(bed, genome, isBED12, mapQual, uncompressedBam);
    }
    else {
        bedtobam_help();
    }
    return 0;
}


void bedtobam_help(void) {
    
    cerr << "\nTool:    bedtools bedtobam (aka bedToBam)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Converts feature records to BAM format." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-mapq\t" << "Set the mappinq quality for the BAM records." << endl;
    cerr                    << "\t\t(INT) Default: 255" << endl << endl;

    cerr << "\t-bed12\t"    << "The BED file is in BED12 format.  The BAM CIGAR" << endl;
    cerr                    << "\t\tstring will reflect BED \"blocks\"." << endl << endl;

    cerr << "\t-ubam\t"     << "Write uncompressed BAM output. Default writes compressed BAM." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  BED files must be at least BED4 to create BAM (needs name field)." << endl << endl;


    // end the program here
    exit(1);
}


void ProcessBed(BedFile *bed, GenomeFile *genome, bool isBED12, int mapQual, bool uncompressedBam) {

    BamWriter *writer = new BamWriter();

    // build a BAM header from the genomeFile
    RefVector refs;
    string    bamHeader;
    map<string, int, std::less<string> > chromToId;
    MakeBamHeader(genome->getGenomeFileName(), refs, bamHeader, chromToId);
        
    // set compression mode
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
    if ( uncompressedBam ) compressionMode = BamWriter::Uncompressed;
    writer->SetCompressionMode(compressionMode);
    // open a BAM and add the reference headers to the BAM file
    writer->Open("stdout", bamHeader, refs);


    // process each BED entry and convert to BAM
    BED bedEntry;
    // open the BED file for reading.
    bed->Open();
    while (bed->GetNextBed(bedEntry)) {
        if (bed->_status == BED_VALID) {
            BamAlignment bamEntry;
            if (bed->bedType >= 4) {
                ConvertBedToBam(bedEntry, bamEntry, chromToId, isBED12, mapQual, bed->_lineNum);
                writer->SaveAlignment(bamEntry);
            }
            else {
                cerr << "Error: BED entry without name found at line: " << bed->_lineNum << ".  Exiting!" << endl;
                exit (1);
            }
        }
    }
    //close up
    bed->Close();
    writer->Close();
}


void ConvertBedToBam(const BED &bed, BamAlignment &bam, map<string, int, std::less<string> > &chromToId,
                     bool isBED12, int mapQual, int lineNum) {

    bam.Name       = bed.name;
    bam.Position   = bed.start;
    bam.Bin        = bedtobam_reg2bin(bed.start, bed.end);

    // hard-code the sequence and qualities.
    int bedLength  = bed.end - bed.start;

    // set dummy seq and qual strings.  the input is BED,
    // so the sequence is inherently the same as it's
    // reference genome.
    // Thanks to James M. Ward for pointing this out.
    bam.QueryBases = "";
    bam.Qualities  = "";

    // chrom and map quality
    bam.RefID      = chromToId[bed.chrom];
    bam.MapQuality = mapQual;

    // set the BAM FLAG
    bam.AlignmentFlag = 0;
    if (bed.strand == "-")
        bam.SetIsReverseStrand(true);

    bam.MatePosition = -1;
    bam.InsertSize   = 0;
    bam.MateRefID    = -1;

    bam.CigarData.clear();

    if (isBED12 == false) {
        CigarOp cOp;
        cOp.Type = 'M';
        cOp.Length = bedLength;
        bam.CigarData.push_back(cOp);
    }
    // we're being told that the input is BED12.
    else{

        // does it smell like BED12?  if so, process it.
        if (bed.fields.size() == 12) {

            // extract the relevant BED fields to convert BED12 to BAM
            // namely: blockCount, blockStarts, blockEnds
            unsigned int blockCount = atoi(bed.fields[9].c_str());

            vector<int> blockSizes, blockStarts;
            Tokenize(bed.fields[10], blockSizes, ',');
            Tokenize(bed.fields[11], blockStarts, ',');

            // make sure this is a well-formed BED12 entry.
            if (blockSizes.size() != blockCount) {
                cerr << "Error: Number of BED blocks does not match blockCount at line: " << lineNum << ".  Exiting!" << endl;
                exit (1);
            }
            else {
                // does the first block start after the bed.start?
                // if so, we need to do some "splicing"
                if (blockStarts[0] > 0) {
                    CigarOp cOp;
                    cOp.Length = blockStarts[0];
                    cOp.Type = 'N';
                    bam.CigarData.push_back(cOp);
                }
                // handle the "middle" blocks
                for (unsigned int i = 0; i < blockCount - 1; ++i) {
                    CigarOp cOp;
                    cOp.Length = blockSizes[i];
                    cOp.Type = 'M';
                    bam.CigarData.push_back(cOp);

                    if (blockStarts[i+1] > (blockStarts[i] + blockSizes[i])) {
                        CigarOp cOp;
                        cOp.Length = (blockStarts[i+1] - (blockStarts[i] + blockSizes[i]));
                        cOp.Type = 'N';
                        bam.CigarData.push_back(cOp);
                    }
                }
                // handle the last block.
                CigarOp cOp;
                cOp.Length = blockSizes[blockCount - 1];
                cOp.Type = 'M';
                bam.CigarData.push_back(cOp);
            }
        }
        // it doesn't smell like BED12.  complain.
        else {
            cerr << "You've indicated that the input file is in BED12 format, yet the relevant fields cannot be found.  Exiting." << endl << endl;
            exit(1);
        }
    }
}


void MakeBamHeader(const string &genomeFile, RefVector &refs, string &header,
                   map<string, int, std::less<string> > &chromToId) {

    // make a genome map of the genome file.
    GenomeFile genome(genomeFile);

    header += "@HD\tVN:1.0\tSO:unsorted\n";
    header += "@PG\tID:BEDTools_bedToBam\tVN:V";
    header += VERSION;
    header += "\n";

    int chromId = 0;
    vector<string> chromList = genome.getChromList();
    // ARQ: 23-May-2012.  No need to sort. Allow genome file to
    // drive the order of the chromosomes in the BAM header.
    // sort(chromList.begin(), chromList.end());

    // create a BAM header (@SQ) entry for each chrom in the BEDTools genome file.
    vector<string>::const_iterator genomeItr  = chromList.begin();
    vector<string>::const_iterator genomeEnd  = chromList.end();
    for (; genomeItr != genomeEnd; ++genomeItr) {
        chromToId[*genomeItr] = chromId;
        chromId++;

        // add to the header text
        int size = genome.getChromSize(*genomeItr);
        string chromLine = "@SQ\tSN:" + *genomeItr + "\tAS:" + genomeFile + "\tLN:" + ToString(size) + "\n";
        header += chromLine;

        // create a chrom entry and add it to
        // the RefVector
        RefData chrom;
        chrom.RefName            = *genomeItr;
        chrom.RefLength          = size;
        refs.push_back(chrom);
    }
}


/* Taken directly from the SAMTools spec
calculate bin given an alignment in [beg,end) (zero-based, half-close, half-open) */
int bedtobam_reg2bin(int beg, int end) {
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
    return 0;
}


