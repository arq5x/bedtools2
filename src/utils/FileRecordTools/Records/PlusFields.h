/*
 * BedPlusInterval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef PLUSFIELDS_H_
#define PLUSFIELDS_H_

#include "QuickString.h"
#include <vector>

using namespace std;

class SingleLineDelimTextFileReader;

class PlusFields {
public:

	PlusFields();
	virtual ~PlusFields() {}
	void setNumOffsetFields(int numOffsetFields) { _numOffsetFields = numOffsetFields; }
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void printFields(QuickString &outBuf) const;

	virtual const QuickString &getField(int fieldNum) const;
	virtual size_t size() const { return _fields.size(); }


protected:
	vector<QuickString> _fields;
	int _numOffsetFields; //could be 3 for BedPlus, but GFF has 8 or 9
};



#endif /* PLUSFIELDS_H_ */
