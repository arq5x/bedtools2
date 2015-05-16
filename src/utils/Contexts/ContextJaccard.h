/*
 * ContextJaccard.h
 *
 *  Created on: Apr 24, 2014
 *      Author: nek3d
 */

#ifndef CONTEXTJACCARD_H_
#define CONTEXTJACCARD_H_

#include "ContextIntersect.h"

class ContextJaccard : public ContextIntersect {
public:
  ContextJaccard();
  ~ContextJaccard();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	virtual bool isValidState();

protected:
	bool handle_s();
	bool handle_S();
};

#endif /* CONTEXTJACCARD_H_ */
