/*
 * SampleFile.cpp
 *
 *  Created on: Nov 18, 2013
 *      Author: nek3d
 */

#include "SampleFile.h"
#include "ContextSample.h"
#include "FileRecordMgr.h"
#include "RecordOutputMgr.h"

static const bool SampleRecordLtFn(const Record *rec1, const Record *rec2) {
	return (*rec1 < *rec2);
}

SampleFile::SampleFile(ContextSample *context)
:	_context(context),
	_inputFile(NULL),
	_outputMgr(NULL),
	_numSamples(0),
	_numCurrSamples(0),
	_currRecordNum(0)
{
	_numSamples = context->getNumOutputRecords();
	if (_numSamples == 0) {
		_numSamples = DEFAULT_NUM_SAMPLES;
	}
}

SampleFile::~SampleFile() {

}

bool SampleFile::takeSample()
{
	//we're only operating on one file, so the idx is zero.
	_inputFile =  _context->getFile(0);
	_samples.resize(_numSamples, NULL);


	//Context object takes care of the seed, either user given or randomly
	//generated, and seeds the call to srand with it, so we don't have to
	//here.
	if (!_context->hasConstantSeed()) {
		_context->getUnspecifiedSeed();
	}

	_outputMgr = new RecordOutputMgr();
	_outputMgr->init(_context);


	while (!_inputFile->eof()) {
		Record *record = _inputFile->getNextRecord();
		if (record == NULL) {
			continue;
		}
		if (!keepRecord(record)) {
			_inputFile->deleteRecord(record);
		}
		_currRecordNum++;
	}

	if (_currRecordNum < _numSamples) {
		//die with error;
		cerr << "\n***** ERROR: Input file has fewer records than the requested number of output records. *****" << endl << endl;
		exit(1);
 	}

	//If the output type is BAM, must sort the output records.
	if (_context->getOutputFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		sort(_samples.begin(), _samples.end(), SampleRecordLtFn);
	}
	// Now output all the kept records, then do cleanup.
	for (size_t i=0; i < _numSamples; i++) {
		_outputMgr->printRecord(_samples[i]);
	}
	delete _outputMgr;
	_inputFile->close();
	return true;
}

bool SampleFile::keepRecord(Record *record)
{
	if (!strandComplies(record)) {
		return false;
	}
	if (_numCurrSamples < _numSamples) {
		_samples[_numCurrSamples] = record;
		_numCurrSamples++;
		return true;
	}


	// We need a random number in the range [0, _currRecordNum].
	// Must combine two consective calls to rand()
    // because RAND_MAX is 2^31 (2147483648), whereas
    // the number of input records could be far larger.
    size_t idx = ((((long) rand()) << 31) | rand()) % _currRecordNum;

    if (idx < _numSamples) {
    	//replace old record at idx with this new one.
    	_inputFile->deleteRecord(_samples[idx]);
    	_samples[idx] = record;
    	return true;
    }
    return false;
}

bool SampleFile::strandComplies(const Record * record) {
	if (!_context->getSameStrand()) {
		return true;
	}
	if (_context->getForwardOnly() && record->getStrandVal() == Record::FORWARD) {
		return true;
	}
	if (_context->getReverseOnly() && record->getStrandVal() == Record::REVERSE) {
		return true;
	}
	return false;
}
