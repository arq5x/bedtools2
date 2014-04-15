/*
 * EmptyRecord.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: nek3d
 */


#include "EmptyRecord.h"
#include "SingleLineDelimTextFileReader.h"

//Technically, there is nothing to do here. Just adding the file
//to keep make happy.

EmptyRecord::EmptyRecord() {

}

EmptyRecord::~EmptyRecord() {

}

bool EmptyRecord::initFromFile(FileReader *fileReader)
{
	SingleLineDelimTextFileReader *sldFileReader = static_cast<SingleLineDelimTextFileReader*>(fileReader);
	return initFromFile(sldFileReader);
}


bool EmptyRecord::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	return true;
}

