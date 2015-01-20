/*
 * ContextSpacing.h
 *
 *  Created on: Jan 15, 2015
 *      Author: Aaron Quinlan
 */

#ifndef CONTEXTSPACING_H_
#define CONTEXTSPACING_H_

#include "ContextBase.h"

class ContextSpacing : public ContextBase {
public:
	ContextSpacing();
	~ContextSpacing();
	//bool isValidState();

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);

};

#endif /* CONTEXTSPACING_H_ */


