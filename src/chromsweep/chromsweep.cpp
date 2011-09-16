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
    
    // prime the results pump.
    _qy_lineNum = 0;
    _db_lineNum = 0;
    
    _hits.reserve(1000);
    _cache.reserve(1000);
    
    _bedA->Open();
    _bedB->Open();
    _qy_status = _bedA->GetNextBed(_curr_qy, _qy_lineNum);
    _db_status = _bedB->GetNextBed(_curr_db, _db_lineNum);
}


/*
    Destructor
*/
ChromSweep::~ChromSweep(void) {
}


bool after(const BED &a, const BED &b) {
    return (a.start >= b.end);
}

void ChromSweep::ScanCache() {
    if (_qy_status != BED_INVALID) {
        vector<BED>::iterator c        = _cache.begin();
        while (c != _cache.end())
        {
            if ((_curr_qy.chrom == c->chrom) && !(after(_curr_qy, *c))) {
                if (overlaps(_curr_qy.start, _curr_qy.end, c->start, c->end) > 0) {
                    _hits.push_back(*c);
                }
                ++c;
            }
            else {
                c = _cache.erase(c);
            }
        }
    }
}

// void ChromSweep::ScanCache(const BED &curr_qy, BedLineStatus qy_status, vector<BED> &db_cache, vector<BED> &hits) {
//     if (qy_status != BED_INVALID) {
//         vector<BED>::iterator c        = db_cache.begin();
//         while (c != db_cache.end())
//         {
//             if ((curr_qy.chrom == c->chrom) && !(after(curr_qy, *c))) {
//                 if (overlaps(curr_qy.start, curr_qy.end, c->start, c->end) > 0) {
//                     hits.push_back(*c);
//                 }
//                 ++c;
//             }
//             else {
//                 c = db_cache.erase(c);
//             }
//         }
//     }
// }


void ChromSweep::ChromCheck() 
{
    if ((_curr_qy.chrom == _curr_db.chrom) || (_db_status == BED_INVALID) || (_qy_status == BED_INVALID)) {
        return;
    }
    
    if (_curr_qy.chrom > _curr_db.chrom) {
        while (!_bedB->Empty() && _curr_db.chrom < _curr_qy.chrom)
        {
            _db_status = _bedB->GetNextBed(_curr_db, _db_lineNum);
        }
        _cache.clear();
    }
    else if (_curr_qy.chrom < _curr_db.chrom) {
        // report hits for the remaining queries on this chrom
        string curr_chrom = _curr_qy.chrom;
        while (!_bedA->Empty() && _curr_qy.chrom == curr_chrom)
        {
            ScanCache();
            _results.push(make_pair(_curr_qy, _hits));
            _qy_status = _bedA->GetNextBed(_curr_qy, _qy_lineNum);
            _hits.clear();
        }
        // now fast forward query to catch up to database
        while (!_bedA->Empty() && _curr_qy.chrom < _curr_db.chrom)
        {
            // hits is empty to reflect the fact that no hits are found in catch-up mode
            _results.push(make_pair(_curr_qy, _hits));
            _qy_status = _bedA->GetNextBed(_curr_qy, _qy_lineNum);
        }
        _cache.clear();
    }
}



void ChromSweep::ReportQuery(const BED &query) {
    _bedA->reportBedTab(query);
}


bool ChromSweep::Next(pair<BED, vector<BED> > &next) {
    if (!_bedA->Empty()) {
        // have we changed chromosomes?
        ChromCheck();
        // scan the database cache for hits
        ScanCache();
        // advance the db until we are ahead of the query. update hits and cache as necessary
        while (!_bedB->Empty() && _curr_qy.chrom == _curr_db.chrom && !(after(_curr_db, _curr_qy)))
        {
            if (overlaps(_curr_qy.start, _curr_qy.end, _curr_db.start, _curr_db.end) > 0) {
                _hits.push_back(_curr_db);
            }
            _cache.push_back(_curr_db);
            _db_status = _bedB->GetNextBed(_curr_db, _db_lineNum);
        }
        // report the hits for this query and reset for the next query
        _results.push(make_pair(_curr_qy, _hits));
        _hits.clear();
        _qy_status = _bedA->GetNextBed(_curr_qy, _qy_lineNum);
    }
    if (!_results.empty()) {
        next = _results.front();
        _results.pop();
        return true;
    }
    else {return false;}
}


// void ChromSweep::Sweep() {
// 
//     int qy_lineNum = 0;
//     int db_lineNum = 0;
// 
//     // current feature from each file
//     BED curr_qy, curr_db;
//     
//     // status of the current lines
//     BedLineStatus qy_status, db_status;
//     vector<BED> db_cache;
//     vector<BED> hits;
//     
//     // open the files; get the first line from each
//     _bedA->Open();
//     _bedB->Open();
//     qy_status = _bedA->GetNextBed(_curr_qy, _qy_lineNum);
//     db_status = _bedB->GetNextBed(_curr_db, _db_lineNum);
//     while (!_bedA->Empty()) {
//         // have we changed chromosomes?
//         ChromCheck(curr_qy, curr_db, qy_status, db_status, qy_lineNum, db_lineNum, db_cache, hits);
//         // scan the database cache for hits
//         //db_cache = ScanCache(curr_qy, qy_status, db_cache, hits);
//         ScanCache(curr_qy, qy_status, db_cache, hits);
//         // advance the db until we are ahead of the query. update hits and cache as necessary
//         while (!_bedB->Empty() && 
//                curr_qy.chrom == curr_db.chrom &&
//                !(after(curr_db, curr_qy)))
//         {
//             if (overlaps(curr_qy.start, curr_qy.end, curr_db.start, curr_db.end) > 0) {
//                 hits.push_back(curr_db);
//             }
//             db_cache.push_back(curr_db);
//             db_status = _bedB->GetNextBed(curr_db, db_lineNum);
//         }
//         // report the hits for this query and reset for the next query
//         //ReportHits(curr_qy, hits);
//         _results.push(make_pair(curr_qy, hits));
//         hits.clear();
//         qy_status = _bedA->GetNextBed(curr_qy, qy_lineNum);
//     }
// }
