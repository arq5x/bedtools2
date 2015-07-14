/*
 * ContextComplement.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#include "ContextComplement.h"

ContextComplement::ContextComplement()
{
	setSortedInput(true);
	setUseMergedIntervals(true);
	setExplicitBedOutput(true);

}

ContextComplement::~ContextComplement()
{

}

bool ContextComplement::isValidState()
{
	if (!hasGenomeFile()) {
		_errorMsg = "\n***** ERROR: no -g genome file provided. *****";
		return false;
	}

	return ContextBase::isValidState();
}
