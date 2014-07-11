/*
 * SingleLineDelimTransferBuffer.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef SINGLELINEDELIMTRANSFERBUFFER_H_
#define SINGLELINEDELIMTRANSFERBUFFER_H_

class SingleLineDelimTransferBuffer {
public:
	SingleLineDelimTransferBuffer(int numFields, char delim='\t');
	~SingleLineDelimTransferBuffer();
	bool initFromInput(const char *inBuffer);
	int getNumFields() const { return _numFields; }
	const char *getField(int fieldNum) const;

protected:
	void clear();

	static const int MAX_FIELD_SIZE = 8192; //Given field can not exceed 8K.
	int _numFields;
	char _delimChar;
	char **_fields;
};


#endif /* SINGLELINEDELIMTRANSFERBUFFER_H_ */
