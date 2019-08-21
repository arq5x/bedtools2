/*
 * spacingFile.cpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nek3d
 */
#include "spacingFile.h"

SpacingFile::SpacingFile(ContextSpacing *context)
: ToolBase(context),
  _prevRec(NULL),
  _currRec(NULL)
{

}

SpacingFile::~SpacingFile()
{
	_inputFile->deleteRecord(_prevRec);
	_prevRec = NULL;
}

bool SpacingFile::init()
{
	//we're only operating on one file, so the idx is zero.
	_inputFile =  _context->getFile(0);
	return true;
}


bool SpacingFile::findNext(RecordKeyVector &hits)
{
	while (!_inputFile->eof()) {

		_currRec = _inputFile->getNextRecord();

		// no more records
		if (_currRec == NULL) {
			continue;
		}
		// first record in file.
		if (_prevRec == NULL) {
			_distance.append(".");
		}
		// the meat of the file
		else {
			// _currRecent and _prevRecious records are on the same chromosome.
			if (_currRec->getChrName() == _prevRec->getChrName())
			{
				// do _currRec and _prevRec overlap?
				if (_currRec->sameChromIntersects(_prevRec, false, false, 1E-9, 1E-9, false, false, false))
					_distance.append("-1");
				else
				{
					CHRPOS distance = _currRec->getStartPos() - _prevRec->getEndPos();
					ostringstream s;
					s << distance;
					_distance.append(s.str());
				}
			}
			// we have changed chromosomes
			else if (_currRec->getChrName() != _prevRec->getChrName())
			{
				_distance.append(".");
			}
		}
		hits.setKey(_currRec);
		return true;
	}
	return false;
}

void SpacingFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	outputMgr->printRecord(hits.getKey(), _distance);
}

void SpacingFile::cleanupHits(RecordKeyVector &hits)
{
	_inputFile->deleteRecord(_prevRec);
	_prevRec = _currRec;
	_distance.clear();
	hits.clearAll();
}
