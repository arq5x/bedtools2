/*
 * SpacingFile.h
 *
 *  Created on: Nov 18, 2013
 *      Author: nek3d
 */

#ifndef SPACINGFILE_H_
#define SPACINGFILE_H_

#include "ContextSpacing.h"
#include "Record.h"
#include <vector>

using namespace std;

class FileRecordMgr;
class Context;
class RecordOutputMgr;

class SpacingFile {
public:
	SpacingFile(ContextSpacing *context);
	~SpacingFile();
	bool getSpacing();

private:
	ContextSpacing *_context;
	FileRecordMgr *_inputFile;
	RecordOutputMgr *_outputMgr;


};

#endif /* SPACINGFILE_H_ */
