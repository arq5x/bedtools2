/*
 * ContextSample.h
 *
 *  Created on: Jan 6, 2014
 *      Author: nek3d
 */

#ifndef CONTEXTSAMPLE_H_
#define CONTEXTSAMPLE_H_

#include "ContextBase.h"

class ContextSample : public ContextBase {
public:
	ContextSample();
	~ContextSample();
	bool isValidState();

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);

	bool getSameStrand() const {return _sameStrand; }
	void setSameStrand(bool val) { _sameStrand = val; }
	bool getForwardOnly() const { return _forwardOnly; }
	bool getReverseOnly() const { return _reverseOnly; }

private:
	virtual bool handle_s();
};


#endif /* CONTEXTSAMPLE_H_ */
