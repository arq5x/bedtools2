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
MultiCovBam::MultiCovBam(const vector<string> &bam_files, const string bed_file, 
                         int minQual, bool properOnly,
                         bool keepDuplicates, bool keepFailedQC)
:
_bam_files(bam_files),
_bed_file(bed_file),
_minQual(minQual),
_properOnly(properOnly),
_keepDuplicates(keepDuplicates),
_keepFailedQC(keepFailedQC)
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
    
    if ( !reader.Open(_bam_files) )
    {
        cerr << "Could not open input BAM files." << endl;
        exit(1);
    }
    else
    {
        // attempt to find index files
        reader.LocateIndexes();

        // if index data available for all BAM files, we can use SetRegion
        if ( reader.HasIndexes() ) {
            BED bed;

            _bed->Open();
            // loop through each BED entry, jump to it, 
            // and collect coverage from each BAM
            while (_bed->GetNextBed(bed))
            {
                if (_bed->_status == BED_VALID)
                {
                    // initialize counts for each file to 0
                    vector<int> counts(_bam_files.size(), 0);
                    // get the BAM refId for this chrom.
                    int refId = reader.GetReferenceID(bed.chrom);
                    // set up a BamRegion to which to attempt to jump
                    BamRegion region(refId, (int)bed.start, refId, (int)bed.end);
                    
                    // everything checks out, just iterate through specified region, counting alignments
                    if ( (refId != -1) && (reader.SetRegion(region)) ) {
                        BamAlignment al;
                        while ( reader.GetNextAlignment(al) )
                        {
                            bool duplicate = al.IsDuplicate();
                            bool failedQC  = al.IsFailedQC();
                            if (_keepDuplicates) duplicate = false;
                            if (_keepFailedQC)    failedQC = false;
                            // map qual must exceed minimum
                            if ((al.MapQuality >= _minQual) && (!duplicate) && (!failedQC)) {
                                // ignore if not properly paired and we actually care.
                                if (_properOnly && !al.IsProperPair())
                                    continue;

                                // lookup the offset of the file name and tabulate 
                                //coverage for the appropriate file
                                counts[bamFileMap[al.Filename]]++;
                            }
                        }
                    }
                    // report the cov at this interval for each file and reset
                    _bed->reportBedTab(bed);
                    ReportCounts(counts);
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
