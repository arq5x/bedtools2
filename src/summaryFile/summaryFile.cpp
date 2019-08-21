/*
 * qcFile.cpp
 *
 *  Created on: Jan 2, 2019
 *      Author: Aaron Quinlan
 */
#include <iomanip>
#include "summaryFile.h"


bool compareByLength(const Interval &a, const Interval &b)
{
    return (a.end - a.start) < b.end - b.start;
}


SummaryFile::SummaryFile(ContextSummary *context)
: ToolBase(context),
  _currRec(NULL),
  _total_length(0),
  _total_intervals(0),
  _genomeFile(context->getGenomeFile()),
  _chromList(_genomeFile->getChromList())
{}

SummaryFile::~SummaryFile()
{}

bool SummaryFile::init()
{
	//we're only operating on one file, so the idx is zero.
	_frm = static_cast<FileRecordMergeMgr *>(upCast(_context)->getFile(0));
	return true;
}


bool SummaryFile::findNext(RecordKeyVector &hits)
{
	while (!_frm->eof()) {

		_currRec = _frm->getNextRecord(&hits);

		// no more records
		if (_currRec == NULL) {
			continue;
		}
		// the meat of the file
		else {
			Interval ivl;
			ivl.start = _currRec->getStartPos();
			ivl.end = _currRec->getEndPos();
			_chromData[_currRec->getChrName()].push_back(ivl);
			_total_intervals++;
			_total_length += ivl.end - ivl.start;
			return true;
		}
	}
	for (auto it=_chromData.begin(); it!=_chromData.end(); ++it)
	{
		sort(it->second.begin(), it->second.end(), compareByLength); 
	}
	return false;  // no more records to find
}

void SummaryFile::giveFinalReport(RecordOutputMgr *outputMgr) 
{

	CHRPOS min_length = LLONG_MAX;
	CHRPOS max_length = 0;
	cout << "chrom\tnum_records\ttotal_bp\tchrom_frac_genome\tfrac_all_ivls\t"
	        "frac_all_bp\tmin\tmax\tmean" 
	     << endl;

	// report the per-chromosome stats
	for (auto chromIt=_chromList.begin(); chromIt!=_chromList.end(); ++chromIt)
	{
		map<string, vector<Interval>>::const_iterator chromDataIt = _chromData.find(*chromIt);

		if (chromDataIt != _chromData.end())
		{
			vector<Interval> chrom_intervals = _chromData[*chromIt];
			CHRPOS num_chrom_intervals = chrom_intervals.size();
			
			// chrom and num_records
		 	cout << *chromIt << "\t" << num_chrom_intervals << "\t";
		 	
		 	// total_bp
		 	CHRPOS chrom_total_bp = 0;
		 	for (auto chromDataIt = chrom_intervals.begin(); chromDataIt != chrom_intervals.end(); ++chromDataIt)
		 	{
		 		chrom_total_bp += chromDataIt->end - chromDataIt->start;
		 	}
		 	cout << chrom_total_bp << "\t";

		 	// chrom_frac_genome
		 	double pct_of_genome = (double) _genomeFile->getChromSize(*chromIt) / 
		 	                      (double) _genomeFile->getGenomeSize();
		 	cout << fixed << std::setprecision(5) << pct_of_genome << "\t";  
		 	             
		 	// frac_all_ivls
		 	double frac_all_ivls = (double) num_chrom_intervals / _total_intervals;
		 	cout << fixed << std::setprecision(3) << frac_all_ivls << "\t";

		 	// frac_all_bp
		 	double frac_all_bp = (double) chrom_total_bp / _total_length;
		 	cout << fixed << std::setprecision(3) << frac_all_bp << "\t";

		 	// min, max
		 	CHRPOS min = chrom_intervals[0].end - chrom_intervals[0].start;
		 	CHRPOS max = chrom_intervals[num_chrom_intervals - 1].end - 
		 	             chrom_intervals[num_chrom_intervals - 1].start;
		 	cout << min << "\t" << max << "\t";

		 	// update overall min and max?
		 	if (min < min_length) min_length = min;
		 	if (max > max_length) max_length = max;

		 	// mean
		 	double mean = (double) chrom_total_bp / (double) num_chrom_intervals;
		 	cout << fixed << std::setprecision(3) << mean << "\t";

			cout << endl;
		}
		else if (chromIt->length() > 0)
		{
			// chrom_frac_genome
		 	double pct_of_genome = (double) _genomeFile->getChromSize(*chromIt) / 
		 	                      (double) _genomeFile->getGenomeSize();
		 	cout << fixed << std::setprecision(2) << pct_of_genome << "\t";  
			// report default if no data available for this chromosome
			cout << *chromIt << "\t" << "0\t0\t"
			     << pct_of_genome << "\t"
			     << "0.0\t0.0\t-1\t-1\t-1"
			     << endl;
		}
	}
	// report the overall stats
	cout << "all\t" << _total_intervals << "\t" << _total_length << "\t"
	     << "1.0" << "\t" << "1.0" << "\t"
	     << min_length << "\t" << max_length << "\t";
	// mean
	double mean = (double) _total_length / (double) _total_intervals;
	cout << fixed << std::setprecision(2) << mean << endl;
}
