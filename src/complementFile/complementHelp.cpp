/*
 * subtractMain.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */


#include "CommonHelp.h"

void complement_help(void) {

    cerr << "\nTool:    bedtools complement (aka complementBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Returns the base pair complement of a feature file." << endl << endl;

    cerr << "Usage:   " << "bedtools complement" << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-L\t"       << "Limit output to solely the chromosomes with records in the input file." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\t     <chromName><TAB><chromSize>" << endl << endl;
    cerr << "\tFor example, Human (hg19):" << endl;
    cerr << "\tchr1\t249250621" << endl;
    cerr << "\tchr2\t243199373" << endl;
    cerr << "\t..." << endl;
    cerr << "\tchr18_gl000207_random\t4262" << endl << endl;

    cerr << "Tip 1. Use samtools faidx to create a genome file from a FASTA: " << endl;
    cerr << "\tOne can the samtools faidx command to index a FASTA file." << endl;
    cerr << "\tThe resulting .fai index is suitable as a genome file, " << endl;
    cerr << "\tas bedtools will only look at the first two, relevant columns" << endl;
    cerr << "\tof the .fai file." << endl << endl;
    cerr << "\tFor example:" << endl;
    cerr << "\tsamtools faidx GRCh38.fa" << endl;
    cerr << "\tbedtools complement -i my.bed -g GRCh38.fa.fai" << endl << endl;

    cerr << "Tip 2. Use UCSC Table Browser to create a genome file: " << endl;
    cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
    cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
    cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \\" << endl;
    cerr << "\t\"select chrom, size from hg19.chromInfo\"  > hg19.genome" << endl << endl;



}
