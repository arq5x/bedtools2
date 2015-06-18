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
	bool isValidState();

protected:

};

#endif /* CONTEXTCOMPLEMENT_H_ */
