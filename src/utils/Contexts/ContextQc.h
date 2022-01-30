/*
 * ContextQc.h
 *
 *  Created on: Jan 2, 2019
 *      Author: Aaron Quinlan
 */

#ifndef CONTEXTQC_H_
#define CONTEXTQC_H_

#include "ContextBase.h"

class ContextQc : public ContextBase {
public:
	ContextQc();
	~ContextQc();
	bool isValidState();

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
};


#endif /* CONTEXTQC_H_ */