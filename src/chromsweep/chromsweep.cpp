/*****************************************************************************
  chromsweep.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "chromsweep.h"
#include <queue>
#include <set>

bool after(const BED &a, const BED &b);
void report_hits(const BED &curr_qy, const vector<BED> &hits);
vector<BED> scan_cache(const BED &curr_qy, BedLineStatus qy_status, const vector<BED> &db_cache, vector<BED> &hits);


/*
    Constructor
*/
ChromSweep::ChromSweep(string bedAFile, string bedBFile, bool anyHit,
                           bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                           float overlapFraction, bool noHit, bool writeCount, bool forceStrand,
                           bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput) {

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _anyHit              = anyHit;
    _noHit               = noHit;
    _writeA              = writeA;
    _writeB              = writeB;
    _writeOverlap        = writeOverlap;
    _writeAllOverlap     = writeAllOverlap;
    _writeCount          = writeCount;
    _overlapFraction     = overlapFraction;
    _forceStrand         = forceStrand;
    _reciprocal          = reciprocal;
    _obeySplits          = obeySplits;
    _bamInput            = bamInput;
    _bamOutput           = bamOutput;

    if (_anyHit || _noHit || _writeCount)
        _printable = false;
    else
        _printable = true;

    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);

    Sweep();
}


/*
    Destructor
*/
ChromSweep::~ChromSweep(void) {
}


bool after(const BED &a, const BED &b) {
    return (a.start >= b.end);
}

void ChromSweep::ScanCache(const BED &curr_qy, BedLineStatus qy_status, vector<BED> &db_cache, vector<BED> &hits) {
    if (qy_status != BED_INVALID) {
        vector<BED>::iterator c        = db_cache.begin();
        while (c != db_cache.end())
        {
            if ((curr_qy.chrom == c->chrom) && !(after(curr_qy, *c))) {
                if (overlaps(curr_qy.start, curr_qy.end, c->start, c->end) > 0) {
                    hits.push_back(*c);
                }
                ++c;
            }
            else {
                c = db_cache.erase(c);
            }
        }
    }
}


void ChromSweep::ChromCheck(BED &curr_qy, BED &curr_db, 
                            BedLineStatus &qy_status, BedLineStatus &db_status,
                            int &qy_lineNum, int &db_lineNum,
                            vector<BED> &db_cache, vector<BED> &hits) 
{
    if ((curr_qy.chrom == curr_db.chrom) || (db_status == BED_INVALID) || (qy_status == BED_INVALID)) {
        return;
    }
    
    if (curr_qy.chrom > curr_db.chrom) {
        while (!_bedB->Empty() && curr_db.chrom < curr_qy.chrom)
        {
            db_status = _bedB->GetNextBed(curr_db, db_lineNum);
        }
        db_cache.clear();
    }
    else if (curr_qy.chrom < curr_db.chrom) {
        // report hits for the remaining queries on this chrom
        BED tmp_curr_qy = curr_qy;
        while (!_bedA->Empty() && tmp_curr_qy.chrom == curr_qy.chrom)
        {
            //db_cache = ScanCache(tmp_curr_qy, qy_status, db_cache, hits);
            ScanCache(tmp_curr_qy, qy_status, db_cache, hits);

            ReportHits(tmp_curr_qy, hits);
            qy_status = _bedA->GetNextBed(tmp_curr_qy, qy_lineNum);
            hits.clear();
        }
        // now fast forward query to catch up to database
        while (!_bedA->Empty() && tmp_curr_qy.chrom < curr_db.chrom)
        {
            // hits is empty to reflect the fact that no hits are found in catch-up mode
            ReportHits(tmp_curr_qy, hits);
            qy_status = _bedA->GetNextBed(tmp_curr_qy, qy_lineNum);
        }
        curr_qy = tmp_curr_qy;
        db_cache.clear();
    }
}


void ChromSweep::ReportHits(const BED &curr_qy, const vector<BED> &hits) {
    _bedA->reportBedTab(curr_qy);
    cout << hits.size() << endl;
}


void ChromSweep::Sweep() {

    int qy_lineNum = 0;
    int db_lineNum = 0;

    // current feature from each file
    BED curr_qy, curr_db;
    
    // status of the current lines
    BedLineStatus qy_status, db_status;
    vector<BED> db_cache;
    vector<BED> hits;
    
    // open the files; get the first line from each
    _bedA->Open();
    _bedB->Open();
    qy_status = _bedA->GetNextBed(curr_qy, qy_lineNum);
    db_status = _bedB->GetNextBed(curr_db, db_lineNum);
    while (!_bedA->Empty()) {
        // have we changed chromosomes?
        ChromCheck(curr_qy, curr_db, qy_status, db_status, qy_lineNum, db_lineNum, db_cache, hits);
        // scan the database cache for hits
        //db_cache = ScanCache(curr_qy, qy_status, db_cache, hits);
        ScanCache(curr_qy, qy_status, db_cache, hits);
        // advance the db until we are ahead of the query. update hits and cache as necessary
        while (!_bedB->Empty() && 
               curr_qy.chrom == curr_db.chrom &&
               !(after(curr_db, curr_qy)))
        {
            if (overlaps(curr_qy.start, curr_qy.end, curr_db.start, curr_db.end) > 0) {
                hits.push_back(curr_db);
            }
            db_cache.push_back(curr_db);
            db_status = _bedB->GetNextBed(curr_db, db_lineNum);
        }
        // report the hits for this query and reset for the next query
        ReportHits(curr_qy, hits);
        hits.clear();
        qy_status = _bedA->GetNextBed(curr_qy, qy_lineNum);
    }
}
