/*****************************************************************************
  mapBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "mapFile.h"
#include "ContextMap.h"
#include "FileRecordMgr.h"
#include "NewChromsweep.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"

const int PRECISION = 21;

FileMap::FileMap(ContextMap *context)
: _context(context),
  _blockMgr(NULL),
  _recordOutputMgr(NULL),
  _colOps(_context->getColOps())
{
  _blockMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
  _recordOutputMgr = new RecordOutputMgr();
  _recordOutputMgr->init(_context);
  _keyListOps.setNullValue(_context->getNullValue());
  _keyListOps.setDelimStr(_context->getDelim());
}

FileMap::~FileMap(void) {
	delete _blockMgr;
	_blockMgr = NULL;
	delete _recordOutputMgr;
	_recordOutputMgr = NULL;
}

bool FileMap::mapFiles()
{
    NewChromSweep sweep(_context);
    if (!sweep.init()) {
      return false;
    }
    RecordKeyList hitSet;
    while (sweep.next(hitSet)) {
    	_outputValues.clear();
    	if (_context->getObeySplits()) {
			RecordKeyList keySet(hitSet.getKey());
			RecordKeyList resultSet(hitSet.getKey());
			_blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
			calculateOutput(resultSet);
			_recordOutputMgr->printRecord(resultSet.getKey(), _outputValues);
    	} else {
			calculateOutput(hitSet);
			_recordOutputMgr->printRecord(hitSet.getKey(), _outputValues);
		}
    }
    return true;
}

void FileMap::calculateOutput(RecordKeyList &hits)
{
	//loop through all requested columns, and for each one, call the method needed
	//for the operation specified.
	_keyListOps.setKeyList(&hits);

	double val = 0.0;
	for (int i=0; i < (int)_colOps.size(); i++) {
		int col = _colOps[i].first;
		KeyListOps::OP_TYPES opCode = _colOps[i].second;

		_keyListOps.setColumn(col);
		switch (opCode) {
		case KeyListOps::SUM:
			val = _keyListOps.getSum();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::MEAN:
			val = _keyListOps.getMean();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::STDDEV:
			val = _keyListOps.getStddev();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::SAMPLE_STDDEV:
			val = _keyListOps.getSampleStddev();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::MEDIAN:
			val = _keyListOps.getMedian();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::MODE:
			_outputValues.append(_keyListOps.getMode());
			break;

		case KeyListOps::ANTIMODE:
			_outputValues.append(_keyListOps.getAntiMode());
			break;

		case KeyListOps::MIN:
			val = _keyListOps.getMin();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::MAX:
			val = _keyListOps.getMax();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::ABSMIN:
			val = _keyListOps.getAbsMin();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::ABSMAX:
			val = _keyListOps.getAbsMax();
			if (isnan(val)) {
				_outputValues.append(_context->getNullValue());
			} else {
				_outputValues.append(val);
			}
			break;

		case KeyListOps::COUNT:
			_outputValues.append(_keyListOps.getCount());
			break;

		case KeyListOps::DISTINCT:
			_outputValues.append(_keyListOps.getDistinct());
			break;

		case KeyListOps::COUNT_DISTINCT:
			_outputValues.append(_keyListOps.getCountDistinct());
			break;

		case KeyListOps::DISTINCT_ONLY:
			_outputValues.append(_keyListOps.getDistinctOnly());
			break;

		case KeyListOps::COLLAPSE:
			_outputValues.append(_keyListOps.getCollapse());
			break;

		case KeyListOps::CONCAT:
			_outputValues.append(_keyListOps.getConcat());
			break;

		case KeyListOps::FREQ_ASC:
			_outputValues.append(_keyListOps.getFreqAsc());
			break;

		case KeyListOps::FREQ_DESC:
			_outputValues.append(_keyListOps.getFreqDesc());
			break;

		case KeyListOps::FIRST:
			_outputValues.append(_keyListOps.getFirst());
			break;

		case KeyListOps::LAST:
			_outputValues.append(_keyListOps.getLast());
			break;

		case KeyListOps::INVALID:
		default:
			// Any unrecognized operation should have been handled already in the context validation.
			// It's thus unnecessary to handle it here, but throw an error to help us know if future
			// refactoring or code changes accidentally bypass the validation phase.
			cerr << "ERROR: Invalid operation given for column " << col << ". Exiting..." << endl;
			break;
		}
		//if this isn't the last column, add a tab.
		if (i < (int)_colOps.size() -1) {
			_outputValues.append('\t');
		}
	}
}
