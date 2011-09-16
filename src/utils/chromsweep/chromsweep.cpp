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

bool after(const BED &a, const BED &b);
void report_hits(const BED &curr_qy, const vector<BED> &hits);
vector<BED> scan_cache(const BED &curr_qy, BedLineStatus qy_status, const vector<BED> &db_cache, vector<BED> &hits);


/*
    // constructor using existing BedFile pointers
*/
ChromSweep::ChromSweep(BedFile *bedA, BedFile *bedB) 
: _bedA(bedA)
, _bedB(bedB)
{
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
    Constructor with filenames
*/
ChromSweep::ChromSweep(string &bedAFile, string &bedBFile) 
{
    // prime the results pump.
    _qy_lineNum = 0;
    _db_lineNum = 0;
    
    _hits.reserve(1000);
    _cache.reserve(1000);
    
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);
    
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


void ChromSweep::ChromCheck() 
{
    // the files are on the chrom
    if ((_curr_qy.chrom == _curr_db.chrom) || (_db_status == BED_INVALID) || (_qy_status == BED_INVALID)) {
        return;
    }
    
    // the query is ahead of the database.
    // fast-forward the database to catch-up.
    if (_curr_qy.chrom > _curr_db.chrom) {
        while (!_bedB->Empty() && _curr_db.chrom < _curr_qy.chrom)
        {
            _db_status = _bedB->GetNextBed(_curr_db, _db_lineNum);
        }
        _cache.clear();
    }
    // the database is ahead of the query.
    // 1. scan the cache for remaining hits on the query's current chrom.
    // 2. fast-forward until we catch up and report 0 hits until we do.
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
        // add the hits for this query to the pump
        _results.push(make_pair(_curr_qy, _hits));
        // reset for the next query
        _hits.clear();
        _curr_qy = _nullBed;
        _qy_status = _bedA->GetNextBed(_curr_qy, _qy_lineNum);
    }
    // report the next set if hits if there are still overlaps in the pump
    if (!_results.empty()) {
        next = _results.front();
        _results.pop();
        return true;
    }
    // otherwise, the party is over.
    else {return false;}
}

