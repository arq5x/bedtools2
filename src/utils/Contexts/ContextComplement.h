/*
 * ContextComplement.h
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#ifndef CONTEXTCOMPLEMENT_H_
#define CONTEXTCOMPLEMENT_H_


#include "ContextBase.h"

class ContextComplement : public ContextBase {
public:
	ContextComplement();
	virtual ~ContextComplement();
  virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	bool isValidState();
  bool getOnlyChromsWithBedRecords() const { return _onlyChromsWithBedRecords; }

private:
  bool handle_limit_chroms();
  bool _onlyChromsWithBedRecords;
protected:

};

#endif /* CONTEXTCOMPLEMENT_H_ */
