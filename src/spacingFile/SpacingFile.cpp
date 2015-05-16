/*
 * SpacingFile.cpp
 *
 *  Created on: Nov 18, 2013
 *      Author: nek3d
 */

#include "SpacingFile.h"
#include "ContextSpacing.h"
#include "FileRecordMgr.h"
#include "RecordOutputMgr.h"


SpacingFile::SpacingFile(ContextSpacing *context)
:	_context(context),
	_inputFile(NULL),
	_outputMgr(NULL)
{}

SpacingFile::~SpacingFile() {

}

bool SpacingFile::getSpacing()
{
	//we're only operating on one file, so the idx is zero.
	_inputFile =  _context->getFile(0);


	_outputMgr = new RecordOutputMgr();
	_outputMgr->init(_context);

	Record *prev = NULL;
	Record *curr = NULL;
	QuickString distance;
	while (!_inputFile->eof()) {

		Record *curr = _inputFile->getNextRecord();
		
		// no more records
		if (curr == NULL) {
			continue;
		}
		// first record in file.
		if (prev == NULL) {
			distance.append(".");
		}
		// the meat of the file
		else {
			// current and previous records are on the same chromosome.
			if (curr->getChrName() == prev->getChrName())
			{
				// do curr and prev overlap?
				if (curr->sameChromIntersects(prev, false, false, 1E-9, false))
					distance.append(-1);
				else
					distance.append(curr->getStartPos() - prev->getEndPos());
			}
			// we have changed chromosomes
			else if (curr->getChrName() != prev->getChrName())
			{
				distance.append(".");
			}
		}
		// report distance between current and previous intervals
		_outputMgr->printRecord(curr, distance);
		// memory cleanup
		_inputFile->deleteRecord(prev);
		distance.clear();

		// current become previous in prep for next record in file.
		prev = curr;
	}

	// cleanup
	delete _outputMgr;
	_inputFile->close();
	return true;
}


