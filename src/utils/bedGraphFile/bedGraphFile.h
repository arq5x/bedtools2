/*****************************************************************************
  bedGraphFile.cpp

  (c) 2010 - Assaf Gordon
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef BEDGRAPHFILE_H
#define BEDGRAPHFILE_H

#include "gzstream.h"
#include "lineFileUtilities.h"
#include "fileType.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include <stdint.h>
#include <cstdio>

using namespace std;

//*************************************************
// Data type tydedef
//*************************************************
#ifndef CHRPOS
typedef int64_t CHRPOS;
#endif

#ifndef DEPTH
typedef uint32_t DEPTH;
#endif

/*
   Structure for regular BedGraph records
 */
template <typename T>
class BEDGRAPH
{
public:
    std::string chrom;
    CHRPOS start;
    CHRPOS end;
    T depth;

public:
    typedef T DEPTH_TYPE;
    // constructors

    // Null
    BEDGRAPH() :
        start(0),
        end(0),
        depth(T())
    {}

    // BEDGraph
    BEDGRAPH(string _chrom, CHRPOS _start, CHRPOS _end, T _depth) :
        chrom(_chrom),
        start(_start),
        end(_end),
        depth(_depth)
    {}
}; // BEDGraph

typedef BEDGRAPH<int32_t> BEDGRAPH_INT;
typedef BEDGRAPH<std::string> BEDGRAPH_STR;
typedef BEDGRAPH<double> BEDGRAPH_FLOAT;

template <typename T>
std::ostream& operator<< (std::ostream& strm, const BEDGRAPH<T>& bg)
{
    strm << bg.chrom << "\t"
        << bg.start << "\t"
        << bg.end << "\t"
        << bg.depth;
    return strm;
}

// enum to flag the state of a given line in a BEDGraph file.
enum BedGraphLineStatus
{
    BEDGRAPH_INVALID = -1,
    BEDGRAPH_HEADER  = 0,
    BEDGRAPH_BLANK   = 1,
    BEDGRAPH_VALID   = 2
};


//************************************************
// BedGraphFile Class methods and elements
//************************************************
class BedGraphFile {

public:

    // Constructor
    BedGraphFile(string &);

    // Destructor
    ~BedGraphFile(void);

    // Open a BEDGraph file for reading (creates an istream pointer)
    void Open(void);

    // Close an opened BED file.
    void Close(void);

    // Get the next BED entry in an opened BED file.
    template <typename T>
    BedGraphLineStatus GetNextBedGraph (BEDGRAPH<T> &bedgraph, int &lineNum)
    {
        // make sure there are still lines to process.
        // if so, tokenize, validate and return the BED entry.
        if (_bedGraphStream->good()) {
            string bedGraphLine;
            vector<string> bedGraphFields;

            // parse the bedStream pointer
            getline(*_bedGraphStream, bedGraphLine);
            if (_bedGraphStream->eof())
                return BEDGRAPH_INVALID;
            if (_bedGraphStream->bad()) {
                cerr << "Error while reading file '" << bedGraphFile << "' : "
                    << strerror(errno) << endl;
                exit(1);
            }
            lineNum++;

            // split into a string vector.
            Tokenize(bedGraphLine,bedGraphFields);
            if (bedGraphLine[bedGraphLine.size()-1] == '\r') {
                bedGraphLine.resize(bedGraphLine.size()-1);
            }

            // load the BED struct as long as it's a valid BED entry.
            return parseLine(bedgraph, bedGraphFields, lineNum);
        }

        // default if file is closed or EOF
        return BEDGRAPH_INVALID;
    }

    // the bedfile with which this instance is associated
    string bedGraphFile;

private:
    // data
    istream *_bedGraphStream;

    template <typename T>
    BedGraphLineStatus parseLine (BEDGRAPH<T> &bg, const vector<string> &lineVector, int &lineNum)
    {
        if (lineVector.size() == 0)
            return BEDGRAPH_BLANK;

        if (lineVector[0].find("track")   != string::npos ||
            lineVector[0].find("browser") != string::npos ||
            lineVector[0].find("#") != string::npos)
            return BEDGRAPH_HEADER;

        if (lineVector.size() != 4)
            return BEDGRAPH_INVALID;

        bg.chrom = lineVector[0];

        stringstream str_start(lineVector[1]);
        if (! (str_start >> bg.start) ) {
            cerr << "Input error, failed to extract start value from '" << lineVector[1]
                << "' (column 2) in " << bedGraphFile << " line " << lineNum << endl;
            exit(1);
        }

        stringstream str_end(lineVector[2]);
        if (! (str_end >> bg.end) ) {
            cerr << "Input error, failed to extract end value from '" << lineVector[2]
                << "' (column 3) in " << bedGraphFile << " line " << lineNum << endl;
            exit(1);
        }

        stringstream str_depth(lineVector[3]);
        if (! (str_depth >> bg.depth) ) {
            cerr << "Input error, failed to extract depth value from '" << lineVector[3]
                << "' (column 4) in " << bedGraphFile << " line " << lineNum << endl;
            exit(1);
        }

        return BEDGRAPH_VALID;
    }
};

#endif /* BEDFILE_H */
