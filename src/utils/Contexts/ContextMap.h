/*
 * ContextMap.h
 *
 *  Created on: Jan 7, 2014
 *      Author: arq5x
 */

#ifndef CONTEXTMAP_H_
#define CONTEXTMAP_H_

#include "ContextIntersect.h"

class ContextMap : public ContextIntersect {
public:
	ContextMap();
	virtual ~ContextMap();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
    virtual bool hasIntersectMethods() const { return true; }
    virtual bool isValidState();


private:

};

#endif /* CONTEXTMAP_H_ */
