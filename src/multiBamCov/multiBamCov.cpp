/*****************************************************************************
  multiBamCov.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "multiBamCov.h"
#include "api/BamMultiReader.h"


/*
    Constructor
*/
MultiCovBam::MultiCovBam(const vector<string> &bam_files, const string bed_file, int minQual, bool properOnly) 
:
_bam_files(bam_files),
_bed_file(bed_file),
_minQual(minQual),
_properOnly(properOnly)
{
	_bed = new BedFile(_bed_file);
    LoadBamFileMap();
}


/*
    Destructor
*/
MultiCovBam::~MultiCovBam(void) 
{}



void MultiCovBam::CollectCoverage()
{
    BamMultiReader reader;
    reader.SetIndexCacheMode(BamIndex::NoIndexCaching);
    
    if ( !reader.Open(_bam_files) )
    {
        cerr << "Could not open input BAM files." << endl; return;
    }
    else
    {
        // attempt to find index files
        reader.LocateIndexes();

        // if index data available for all BAM files, we can use SetRegion
        if ( reader.HasIndexes() ) {
            BED bed, nullBed;
            int lineNum = 0;
            BedLineStatus bedStatus;

            _bed->Open();
            // loop through each BED entry, jump to it, 
            // and collect coverage from each BAM
            while ((bedStatus = _bed->GetNextBed(bed, lineNum)) != BED_INVALID)
            {
                if (bedStatus == BED_VALID)
                {
                    reader.SetRegion(reader.GetReferenceID(bed.chrom),
                                           (int) bed.start,
                                           reader.GetReferenceID(bed.chrom),
                                           (int) bed.end);
                                           
                    // everything checks out, just iterate through specified region, counting alignments
                    vector<int> counts(_bam_files.size());
                    BamAlignment al;
                    while ( reader.GetNextAlignment(al) )
                    {
                        // map qual must exceed minimum
                        if (al.MapQuality >= _minQual) {
                            // ignore if not properly paired and we actually care.
                            if (_properOnly && !al.IsProperPair())
                                continue;

                            // lookup the offset of the file name and tabulate 
                            //coverage for the appropriate file
                            counts[bamFileMap[al.Filename]]++;
                        }
                    }
                    // report the cov at this interval for each file and reset
                    _bed->reportBedTab(bed);
                    ReportCounts(counts);
                    bed = nullBed;
                }
            }
            _bed->Close();
        }
        else {
            cerr << "Could not find indexes." << endl;
            reader.Close();
            exit(1);
        }
    }
}


void MultiCovBam::LoadBamFileMap(void) 
{
    for (size_t i = 0; i < _bam_files.size(); ++i)
    {
        bamFileMap[_bam_files[i]] = i;
    }
}

void MultiCovBam::ReportCounts(const vector<int> &counts) 
{
    for (size_t i = 0; i < counts.size(); ++i)
    {
        if (i < counts.size() - 1)
            cout << counts[i] << "\t";
        else
            cout << counts[i];
    }
    cout << endl;
}
