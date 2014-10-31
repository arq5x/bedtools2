/*
 * newClosestMain.cpp
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

using namespace std;

#include "closestFile.h"
#include "ContextClosest.h"

// define our program name
#define PROGRAM_NAME "bedtools closest"

void closest_help(void);

int closest_main(int argc, char* argv[]) {

    ContextClosest *context = new ContextClosest();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
        if (!context->getErrorMsg().empty()) {
            cerr << context->getErrorMsg() << endl;
        }
        closest_help();
        delete context;
        return 0;
    }
    ClosestFile *closestFile = new ClosestFile(context);

    bool retVal = closestFile->getClosest();
    delete closestFile;
    delete context;
    return retVal ? 0 : 1;
}

void closest_help(void) {

    cerr << "\nTool:    bedtools closest (aka closestBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: For each feature in A, finds the closest " << endl;
    cerr << "\t feature (upstream or downstream) in B." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"            << "Req. same strandedness.  That is, find the closest feature in" << endl;
    cerr                        << "\t\tB that overlaps A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Req. opposite strandedness.  That is, find the closest feature" << endl;
    cerr                        << "\t\tin B that overlaps A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-d\t"            << "In addition to the closest feature in B, " << endl;
    cerr                        << "\t\treport its distance to A as an extra column." << endl;
    cerr                        << "\t\t- The reported distance for overlapping features will be 0." << endl << endl;

    cerr << "\t-D\t"            << "Like -d, report the closest feature in B, and its distance to A" << endl;
    cerr                        << "\t\tas an extra column. Unlike -d, use negative distances to report" << endl;
    cerr                        << "\t\tupstream features." << endl;
    cerr                        << "\t\tThe options for defining which orientation is \"upstream\" are:" << endl;
    cerr                        << "\t\t- \"ref\"   Report distance with respect to the reference genome. " << endl;
    cerr                        << "\t\t            B features with a lower (start, stop) are upstream" << endl;
    cerr                        << "\t\t- \"a\"     Report distance with respect to A." << endl;
    cerr                        << "\t\t            When A is on the - strand, \"upstream\" means B has a" << endl;
    cerr                        << "\t\t            higher (start,stop)." << endl;
    cerr                        << "\t\t- \"b\"     Report distance with respect to B." << endl;
    cerr                        << "\t\t            When B is on the - strand, \"upstream\" means A has a" << endl;
    cerr                        << "\t\t            higher (start,stop)." << endl << endl;

    cerr << "\t-io\t"           << "Ignore features in B that overlap A.  That is, we want close," << endl;
    cerr                        << "\t\tyet not touching features only." << endl << endl;

    cerr << "\t-iu\t"           << "Ignore features in B that are upstream of features in A." << endl;
    cerr                        << "\t\tThis option requires -D and follows its orientation" << endl;
    cerr                        << "\t\trules for determining what is \"upstream\"." << endl;
    cerr << "\t-id\t"           << "Ignore features in B that are downstream of features in A." << endl;
    cerr                        << "\t\tThis option requires -D and follows its orientation" << endl;
    cerr                        << "\t\trules for determining what is \"downstream\"." << endl;

    cerr << "\t-t\t"            << "How ties for closest feature are handled.  This occurs when two" << endl;
    cerr                        << "\t\tfeatures in B have exactly the same \"closeness\" with A." << endl;
    cerr                        << "\t\tBy default, all such features in B are reported." << endl;
    cerr                        << "\t\tHere are all the options:" << endl;
    cerr                        << "\t\t- \"all\"    Report all ties (default)." << endl;
    cerr                        << "\t\t- \"first\"  Report the first tie that occurred in the B file." << endl;
    cerr                        << "\t\t- \"last\"   Report the last tie that occurred in the B file." << endl << endl;

    cerr << "\t-mdb\t"          << "How multiple databases are resolved." << endl;
    cerr                        << "\t\t- \"each\"    Report closest records for each database (default)." << endl;
    cerr                        << "\t\t- \"all\"  Report closest records among all databases." << endl << endl;


    cerr << "\t-N\t"            << "Require that the query and the closest hit have different names." << endl;
    cerr                        << "\t\tFor BED, the 4th column is compared." << endl << endl;

    cerr << "\t-header\t"       << "Print the header from the A file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\tReports \"none\" for chrom and \"-1\" for all other fields when a feature" << endl;
    cerr << "\tis not found in B on the same chromosome as the feature in A." << endl;
    cerr << "\tE.g. none\t-1\t-1" << endl << endl;

    // end the program here
    exit(1);
}
