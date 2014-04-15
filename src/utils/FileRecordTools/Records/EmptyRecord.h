/*
 * EmptyRecord.h
 *
 *  Created on: Apr 14, 2014
 *      Author: nek3d
 */

#ifndef EMPTYRECORD_H_
#define EMPTYRECORD_H_

#include "Record.h"

class SingleLineDelimTextFileReader;

class EmptyRecord : public Record {
public:
	//This is really just a place holder that doesn't have to do much.
	//But in order for the rest of the code to not have to contain special
	//ugly if statements for empty records, it at least has to instantiable,
	//so any purely virtual functions from the base class will have to be overriden.
	EmptyRecord();
	~EmptyRecord();

	bool initFromFile(FileReader *);
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	int getNumFields() const { return 0; }
};


#endif /* EMPTYRECORD_H_ */
