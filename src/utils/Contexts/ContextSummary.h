/*
 * ContextSummary.h
 *
 *  Created on: Jan 2, 2019
 *      Author: Aaron Quinlan
 */

#ifndef CONTEXTSUMMARY_H_
#define CONTEXTSUMMARY_H_

#include "ContextBase.h"

class ContextSummary : public ContextBase {
public:
	ContextSummary();
	~ContextSummary();
	bool isValidState();

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
};


#endif /* CONTEXTSUMMARY_H_ */