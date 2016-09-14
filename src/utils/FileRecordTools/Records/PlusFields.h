/*
 * BedPlusInterval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef PLUSFIELDS_H_
#define PLUSFIELDS_H_

#include "string.h"
#include <string>
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
	virtual void printFields(string &outBuf) const;

	virtual const string &getField(int fieldNum) const;
	virtual size_t size() const { return _fields.size(); }


protected:
	vector<string> _fields;
	int _numOffsetFields; //could be 3 for BedPlus, but GFF has 8 or 9
};



#endif /* PLUSFIELDS_H_ */
