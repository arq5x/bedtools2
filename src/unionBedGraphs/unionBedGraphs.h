/*****************************************************************************
  unionBedGraphs.h

  (c) 2010 - Assaf Gordon, CSHL
           - Aaron Quinlan, UVA
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef UNIONBEDGRAPHS_H
#define UNIONBEDGRAPHS_H

#include <vector>
#include <string>
#include "bedGraphFile.h"
#include "genomeFile.h"
#include "Point.h"

class UnionBedGraphs
{
private:
    typedef BEDGRAPH_STR BEDGRAPH_TYPE;

    vector<string>  filenames;
    vector<string>  titles;

    vector<BedGraphFile*>               bedgraph_files;
    vector<BEDGRAPH_TYPE::DEPTH_TYPE>   current_depth;
    vector<BEDGRAPH_TYPE>               current_bedgraph_item;

    std::ostream    &output;

    POINTWITHDEPTH_PQUEUE queue;
    std::string              current_chrom;
    int                      current_non_zero_inputs;
    bool                     print_empty_regions;

    GenomeFile* genome_sizes;

    std::string no_coverage_value;

public:
    UnionBedGraphs(std::ostream& _output,
            const vector<string>& _filenames,
            const vector<string>& _titles,
            bool _print_empty_regions,
            const std::string& _genomeFileName,
            const std::string& _no_coverage_value);

    virtual ~UnionBedGraphs();

    // Combines all bedgraph files
    void Union();

    // Print the header line: chrom/start/end + name of each bedgraph file.
    void PrintHeader();


private:

    // Open all BedGraph files, initialize "current_XXX" vectors
    void OpenBedgraphFiles();

    // Close the BedGraph files.
    void CloseBedgraphFiles();

    /*
       Add an interval from BedGraph file 'index' into the queue.
       will only be added if it belongs to the current chromosome.

       If the interval was added (=consumed), the next interval will be read from the file
       using 'LoadNextBedgraphItem'
     */
    void AddInterval(int index);

    /*
       Loads the next interval from BedGraph file 'index'.
       Stores it in 'current_bedgraph_item' vector.
     */
    void LoadNextBedgraphItem(int index);

    /*
       Scans the 'current_bedgraph_item' vector,
       find the 'first' chromosome to use (different BedGraph files can start with different chromosomes).
     */
    std::string DetermineNextChrom();

    /*
       Returns 'true' if ALL intervals from ALL BedGraph files were used
    */
    bool        AllFilesDone();

    /*
       Extract the next coordinate from the queue, and updates the current coverage information.
       If multiple interval share the same coordinate values, all of them are handled.
       If an END coordinate is consumed, the next interval (from the corresponding file) is read.
     */
    CHRPOS ConsumeNextCoordinate();

    /*
       Updates the coverage information based on the given item.
       Item can be a START coordinate or an END coordiante.
     */
    void UpdateInformation(const PointWithDepth &item);

    /*
       prints chrom/start/end and the current depth coverage values of all the files.
     */
    void PrintCoverage(CHRPOS start, CHRPOS end);

    /*
       prints chrom/start/end and the ZERO depth coverage values of all the files.
     */
    void PrintEmptyCoverage(CHRPOS start, CHRPOS end);

    void DebugPrintQueue();
};


#endif
