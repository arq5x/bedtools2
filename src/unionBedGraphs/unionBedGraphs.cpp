/*****************************************************************************
  unionBedGraphs.cpp

  (c) 2010 - Assaf Gordon, CSHL
           - Aaron Quinlan, UVA
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "bedGraphFile.h"
#include "unionBedGraphs.h"

using namespace std;


UnionBedGraphs::UnionBedGraphs(std::ostream& _output,
                            const vector<string>& _filenames,
                            const vector<string>& _titles,
                            bool _print_empty_regions,
                            const std::string& _genome_size_filename,
                            const std::string& _no_coverage_value   ) :
    filenames(_filenames),
    titles(_titles),
    output(_output),
    current_non_zero_inputs(0),
    print_empty_regions(_print_empty_regions),
    genome_sizes(NULL),
    no_coverage_value(_no_coverage_value)
{
    if (print_empty_regions) {
        assert(!_genome_size_filename.empty());

        genome_sizes = new GenomeFile(_genome_size_filename);
    }
}


UnionBedGraphs::~UnionBedGraphs() {
    CloseBedgraphFiles();
    if (genome_sizes) {
        delete genome_sizes;
        genome_sizes = NULL ;
    }
}


void UnionBedGraphs::Union() {
    OpenBedgraphFiles();

    // Add the first interval from each file
    for(size_t i=0;i<bedgraph_files.size();++i)
        LoadNextBedgraphItem(i);

    // Chromosome loop - once per chromosome
    do {
        // Find the first chromosome to use
        current_chrom = DetermineNextChrom();

        // Populate the queue with initial values from all files
        // (if they belong to the correct chromosome)
        for(size_t i=0;i<bedgraph_files.size();++i)
            AddInterval(i);

        CHRPOS current_start = ConsumeNextCoordinate();

        // User wanted empty regions, and the first coordinate is not 0 - print a dummy empty coverage
        if (print_empty_regions && current_start > 0)
            PrintEmptyCoverage(0,current_start);

        // Intervals loop - until all intervals (of current chromosome) from all files are used.
        do {
            CHRPOS current_end = queue.top().coord;
            PrintCoverage(current_start, current_end);
            current_start = ConsumeNextCoordinate();
        } while (!queue.empty());

        // User wanted empty regions, and the last coordinate is not the last coordinate of the chromosome
            // print a dummy empty coverage
        if (print_empty_regions) {
            CHRPOS chrom_size = genome_sizes->getChromSize(current_chrom);
            if (current_start < chrom_size)
                PrintEmptyCoverage(current_start, chrom_size);
        }

    } while (!AllFilesDone());
}


CHRPOS UnionBedGraphs::ConsumeNextCoordinate() {
    assert(!queue.empty());

    CHRPOS new_position = queue.top().coord;
    do {
        PointWithDepth item = queue.top();
        UpdateInformation(item);
        queue.pop();
    } while (!queue.empty() && queue.top().coord == new_position);

    return new_position;
}


void UnionBedGraphs::UpdateInformation(const PointWithDepth &item) {
    // Update the depth coverage for this file

    // Which coordinate is it - start or end?
    switch (item.coord_type)
    {
    case START:
        current_depth[item.source_index] = item.depth;
        current_non_zero_inputs++;
        break;
    case END:
        //Read the next interval from this file
        AddInterval(item.source_index);
        current_depth[item.source_index] = no_coverage_value;
        current_non_zero_inputs--;
        break;
    default:
        assert(0);
    }
}


void UnionBedGraphs::PrintHeader() {
    output << "chrom\tstart\tend" ;
    for (size_t i=0;i<titles.size();++i)
        output << "\t" <<titles[i];
    output << endl;
}


void UnionBedGraphs::PrintCoverage(CHRPOS start, CHRPOS end) {
    if ( current_non_zero_inputs == 0 && ! print_empty_regions )
        return ;

    output << current_chrom << "\t"
        << start << "\t"
        << end;

    for (size_t i=0;i<current_depth.size();++i)
        output << "\t" << current_depth[i] ;

    output << endl;
}


void UnionBedGraphs::PrintEmptyCoverage(CHRPOS start, CHRPOS end) {
    output << current_chrom << "\t"
        << start << "\t"
        << end;

    for (size_t i=0;i<current_depth.size();++i)
        output << "\t" << no_coverage_value ;

    output << endl;
}


void UnionBedGraphs::LoadNextBedgraphItem(int index) {
    assert(static_cast<unsigned int>(index) < bedgraph_files.size());

    current_bedgraph_item[index].chrom="";

    BedGraphFile *file = bedgraph_files[index];
    BEDGRAPH_STR bg;
    int lineNum = 0;
    BedGraphLineStatus status;

    while ( (status = file->GetNextBedGraph(bg, lineNum)) != BEDGRAPH_INVALID )  {
        if (status != BEDGRAPH_VALID)
            continue;
        current_bedgraph_item[index] = bg;
        break;
    }
}


bool UnionBedGraphs::AllFilesDone() {
    for (size_t i=0;i<current_bedgraph_item.size();++i)
        if (!current_bedgraph_item[i].chrom.empty())
            return false;
    return true;
}


string UnionBedGraphs::DetermineNextChrom() {
    string next_chrom;
    for (size_t i=0;i<current_bedgraph_item.size();++i) {
        if (current_bedgraph_item[i].chrom.empty())
            continue;

        if (next_chrom.empty())
            next_chrom = current_bedgraph_item[i].chrom;
        else
            if (current_bedgraph_item[i].chrom < next_chrom)
                next_chrom = current_bedgraph_item[i].chrom ;
    }
    return next_chrom;
}


void UnionBedGraphs::AddInterval(int index) {
    assert(static_cast<unsigned int>(index) < bedgraph_files.size());

    //This file has no more intervals
    if (current_bedgraph_item[index].chrom.empty())
        return ;

    //If the next interval belongs to a different chrom, don't add it
    if (current_bedgraph_item[index].chrom!=current_chrom)
        return ;

    const BEDGRAPH_STR &bg(current_bedgraph_item[index]);

    PointWithDepth start_item(index, START, bg.start, bg.depth);
    PointWithDepth end_item(index, END, bg.end, bg.depth);

    queue.push(start_item);
    queue.push(end_item);

    LoadNextBedgraphItem(index);
}


void UnionBedGraphs::OpenBedgraphFiles() {
    for (size_t i=0;i<filenames.size();++i) {
        BedGraphFile *file = new BedGraphFile(filenames[i]);
        file->Open();
        bedgraph_files.push_back(file);

        current_depth.push_back(no_coverage_value);
    }
    current_bedgraph_item.resize(filenames.size());
}


void UnionBedGraphs::CloseBedgraphFiles() {
    for (size_t i=0;i<bedgraph_files.size();++i) {
        BedGraphFile *file = bedgraph_files[i];
        delete file;
        bedgraph_files[i] = NULL ;
    }
    bedgraph_files.clear();
}
