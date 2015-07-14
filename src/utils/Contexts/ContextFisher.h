/*
 * ContextFisher.h
 *
 *  Created on: Apr 24, 2014
 *      Author: nek3d
 */

#ifndef CONTEXTFISHER_H_
#define CONTEXTFISHER_H_

#include "ContextJaccard.h"
#include "GenomeFile.h"

class ContextFisher : public ContextJaccard {
public:
  ContextFisher();
  ~ContextFisher();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	virtual bool isValidState();
	string getExcludeFile() { return _excludeFile; }
	void setExcludeFile(string excludeFile) { _excludeFile = excludeFile; }


protected:
	string _excludeFile;
	bool handle_exclude();
	bool handle_m();


};

#endif /* CONTEXTFISHER_H_ */
