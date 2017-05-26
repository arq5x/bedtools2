/*
 * ContextMerge.h
 *
 *  Created on: Mar 26, 2014
 *      Author: nek3d
 */

#ifndef CONTEXTMERGE_H_
#define CONTEXTMERGE_H_
#include <unistd.h>
#include "ContextBase.h"
#include "FileRecordMergeMgr.h"

class ContextMerge: public ContextBase {
public:
	ContextMerge();
	~ContextMerge();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	virtual bool isValidState();

protected:
	bool handle_d();
	bool handle_n();
	bool handle_nms();
	bool handle_scores();
	bool handle_s();
	bool handle_S();
};


#endif /* CONTEXTMERGE_H_ */
