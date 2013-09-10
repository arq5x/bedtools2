/*****************************************************************************
  multiIntersectBed.h

  (c) 2010 - Aaron Quinlan, UVA
           - Assaf Gordon, CSHL
  Quinlan Laboratory
  Department of Public Health Sciences
  Center for Public Health Genomics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef MULTIINTERSECTBED_H
#define MULTIINTERSECTBED_H

#include <vector>
#include <string>
#include "bedFile.h"
#include "GenomeFile.h"
#include "Point.h"

class MultiIntersectBed
{
private:

    vector<string>  filenames;
    vector<string>  titles;

    vector<BedFile*>   input_files;
    vector<int>        current_depth;
    vector<BED>        current_item;

    std::ostream    &output;

    POINT_PQUEUE queue;
    std::string              current_chrom;
    map<int, bool>           files_with_coverage;
    int                      current_non_zero_inputs;
    bool                     print_empty_regions;
    bool                     haveTitles;
    
    GenomeFile* genome_sizes;

    std::string no_coverage_value;

public:
    MultiIntersectBed(std::ostream& _output,
            const vector<string>& _filenames,
            const vector<string>& _titles,
            bool _print_empty_regions,
            const std::string& _genomeFileName,
            const std::string& _no_coverage_value);

    virtual ~MultiIntersectBed();

    // Combines all interval files
    void MultiIntersect();
    
    // 
    void Cluster();

    // Print the header line: chrom/start/end + name of each bedgraph file.
    void PrintHeader();


private:

    // Open all input files, initialize "current_XXX" vectors
    void OpenFiles();

    // Close the input files.
    void CloseFiles();

    /*
       Add an interval from BedGraph file 'index' into the queue.
       will only be added if it belongs to the current chromosome.

       If the interval was added (=consumed), the next interval will be read from the file
       using 'LoadNextItem'
     */
    void AddInterval(int index);

    /*
       Loads the next interval from Bed file 'index'.
       Stores it in 'current_bed_item' vector.
     */
    void LoadNextItem(int index);

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
