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
	string getExcludeFile() { return _excludeFile; }
	void setExcludeFile(string excludeFile) { _excludeFile = excludeFile; }


private:
	bool handle_s();
	bool handle_S();
	bool handle_exclude();
	string _excludeFile;


};

#endif /* CONTEXTFISHER_H_ */
