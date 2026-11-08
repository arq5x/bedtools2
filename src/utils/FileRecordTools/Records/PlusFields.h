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
	PlusFields(const PlusFields &other);
	PlusFields &operator=(const PlusFields &other);
	virtual ~PlusFields();
	void setNumOffsetFields(int numOffsetFields) { _numOffsetFields = numOffsetFields; }
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void printFields(string &outBuf) const;

	virtual const string &getField(int fieldNum) const;
	virtual size_t size() const { return _numFields; }


protected:
	string _text; //the extra columns, tab-joined, exactly as in the line
	size_t _numFields;
	vector<int> _lens; //field lengths, only when the file delimiter is not a tab
	mutable vector<string> *_split; //per-column strings, built on first getField
	int _numOffsetFields; //could be 3 for BedPlus, but GFF has 8 or 9
};



#endif /* PLUSFIELDS_H_ */
