#include "CommonHelp.h"

void groupby_help(void) {

    cerr << "\nTool:    bedtools groupby " << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Summarizes a dataset column based upon" << endl;
    cerr << "\t common column groupings. Akin to the SQL \"group by\" command." << endl << endl;

    cerr << "Usage:\t " << "bedtools groupby" << " -g [group_column(s)] -c [op_column(s)] -o [ops] " << endl;
    cerr << "\t "     << "cat [FILE] | " << "bedtools groupby" << " -g [group_column(s)] -c [op_column(s)] -o [ops] " << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-i\t\t"        << "Input file. Assumes \"stdin\" if omitted." << endl << endl;

    cerr << "\t-g -grp\t\t"      << "Specify the columns (1-based) for the grouping." << endl;
    cerr                         << "\t\t\tThe columns must be comma separated." << endl;
    cerr                         << "\t\t\t- Default: 1,2,3" << endl << endl;

    cerr << "\t-c -opCols\t"     << "Specify the column (1-based) that should be summarized." << endl;
    cerr                         << "\t\t\t- Required." << endl << endl;

    cerr << "\t-o -ops\t\t"      << "Specify the operation that should be applied to opCol." << endl;
    cerr                         << "\t\t\tValid operations:" << endl;
    cerr                         << "\t\t\t    sum, count, count_distinct, min, max," << endl;
    cerr                         << "\t\t\t    mean, median, mode, antimode," << endl;
    cerr                         << "\t\t\t    stdev, sstdev (sample standard dev.)," << endl;
    cerr                         << "\t\t\t    collapse (i.e., print a comma separated list (duplicates allowed)), " << endl;
    cerr                         << "\t\t\t    distinct (i.e., print a comma separated list (NO duplicates allowed)), " << endl;
    cerr                         << "\t\t\t    distinct_sort_num (as distinct, but sorted numerically, ascending), " << endl;
    cerr                         << "\t\t\t    distinct_sort_num_desc (as distinct, but sorted numerically, descending), " << endl;
    cerr                         << "\t\t\t    concat   (i.e., merge values into a single, non-delimited string), " << endl;
    cerr                         << "\t\t\t    freqdesc (i.e., print desc. list of values:freq)" << endl;
    cerr                         << "\t\t\t    freqasc (i.e., print asc. list of values:freq)" << endl;
    cerr                         << "\t\t\t    first (i.e., print first value)" << endl;
    cerr                         << "\t\t\t    last (i.e., print last value)" << endl;

    cerr                         << "\t\t\t- Default: sum" << endl << endl;

    cerr						<< "\t\tIf there is only column, but multiple operations, all operations will be" << endl;
    cerr						<< "\t\tapplied on that column. Likewise, if there is only one operation, but" << endl;
    cerr						<< "\t\tmultiple columns, that operation will be applied to all columns." << endl;
    cerr						<< "\t\tOtherwise, the number of columns must match the the number of operations," << endl;
    cerr						<< "\t\tand will be applied in respective order." << endl;
    cerr						<< "\t\tE.g., \"-c 5,4,6 -o sum,mean,count\" will give the sum of column 5," << endl;
    cerr						<< "\t\tthe mean of column 4, and the count of column 6." << endl;
    cerr						<< "\t\tThe order of output columns will match the ordering given in the command." << endl << endl<<endl;

    cerr << "\t-full\t\t"   << "Print all columns from input file.  The first line in the group is used." << endl;
    cerr            << "\t\t\tDefault: print only grouped columns." << endl << endl;

    cerr << "\t-inheader\t" << "Input file has a header line - the first line will be ignored." << endl << endl ;

    cerr << "\t-outheader\t"    << "Print header line in the output, detailing the column names. " << endl;
    cerr            << "\t\t\tIf the input file has headers (-inheader), the output file" << endl;
    cerr            << "\t\t\twill use the input's column names." << endl;
    cerr            << "\t\t\tIf the input file has no headers, the output file" << endl;
    cerr            << "\t\t\twill use \"col_1\", \"col_2\", etc. as the column names." << endl << endl;

    cerr << "\t-header\t\t" << "same as '-inheader -outheader'" << endl << endl;

    cerr << "\t-ignorecase\t"   << "Group values regardless of upper/lower case." << endl << endl;

    cerr << "\t-prec\t"   << "Sets the decimal precision for output (Default: 5)" << endl << endl;
    cerr << "\t-delim\t"                 << "Specify a custom delimiter for the collapse operations." << endl;
    cerr                                 << "\t\t- Example: -delim \"|\"" << endl;
    cerr                                 << "\t\t- Default: \",\"." << endl << endl;

    cerr << "Examples: " << endl;
    cerr << "\t$ cat ex1.out" << endl;
    cerr << "\tchr1 10  20  A   chr1    15  25  B.1 1000    ATAT" << endl;
    cerr << "\tchr1 10  20  A   chr1    25  35  B.2 10000   CGCG" << endl << endl;
    cerr << "\t$ groupBy -i ex1.out -g 1,2,3,4 -c 9 -o sum" << endl;
    cerr << "\tchr1 10  20  A   11000" << endl << endl;
    cerr << "\t$ groupBy -i ex1.out -grp 1,2,3,4 -opCols 9,9 -ops sum,max" << endl;
    cerr << "\tchr1 10  20  A   11000   10000" << endl << endl;
    cerr << "\t$ groupBy -i ex1.out -g 1,2,3,4 -c 8,9 -o collapse,mean" << endl;
    cerr << "\tchr1 10  20  A   B.1,B.2,    5500" << endl << endl;
    cerr << "\t$ cat ex1.out | groupBy -g 1,2,3,4 -c 8,9 -o collapse,mean" << endl;
    cerr << "\tchr1 10  20  A   B.1,B.2,    5500" << endl << endl;
    cerr << "\t$ cat ex1.out | groupBy -g 1,2,3,4 -c 10 -o concat" << endl;
    cerr << "\tchr1 10  20  A   ATATCGCG" << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  The input file/stream should be sorted/grouped by the -grp. columns" << endl;
    cerr << "\t(2)  If -i is unspecified, input is assumed to come from stdin." << endl << endl;


    // end the program here
    exit(1);

}
