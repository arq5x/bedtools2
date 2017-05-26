/*
 * mapHelp.cpp
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#include "CommonHelp.h"
#include "KeyListOps.h"

void map_help(void) {

    cerr << "\nTool:    bedtools map (aka mapBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Apply a function to a column from B intervals that overlap A." << endl << endl;

    cerr << "Usage:   " << "bedtools map" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;


    KeyListOpsHelp();

    IntersectCommonHelp();
    allToolsCommonHelp();

    cerr << "Notes: " << endl;
    cerr << "\t(1) Both input files must be sorted by chrom, then start." << endl << endl;
}



