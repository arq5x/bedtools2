/*
 * ContextFisher.h
 *
 *  Created on: Apr 24, 2014
 *      Author: nek3d
 */

#ifndef CONTEXTFISHER_H_
#define CONTEXTFISHER_H_

#include "ContextIntersect.h"
#include "GenomeFile.h"

class ContextFisher : public ContextIntersect {
public:
  ContextFisher();
  ~ContextFisher();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	virtual bool isValidState();

private:
	bool handle_s();
	bool handle_S();
};

#endif /* CONTEXTFISHER_H_ */
