/*
 * BamFileReader.h
 *
 *  Created on: Dec 4, 2012
 *      Author: nek3d
 */

#ifndef BAMFILEREADER_H_
#define BAMFILEREADER_H_

#include "FileReader.h"
#include "string.h"
#include "api/BamReader.h"
#include "api/BamAux.h"

class BamFileReader : public FileReader {
public:
	BamFileReader();
	virtual ~BamFileReader();
	virtual bool open(); //open the file
	virtual bool isOpen() const;
	virtual void close();
	virtual bool eof() const {
		return _eof;
	}

	//setUseTags will tell the BamReader to give us all
	//the extra tag information in a BAM record. By default,
	//this is set to false, so not using them, which reduces
	//the run time of reading a BAM file by more than half.
	virtual void setUseTags(bool flag) { _useTags = flag; }
	void setBamReader(BamTools::BamReader *bamReader) { _bamReader = bamReader; }
	virtual bool readEntry();

	virtual bool hasHeader() const { return _bamReader->IsOpen(); } //any open Bam file automatically has a header
	virtual const string &getHeader() const { return _bamHeader; }
	const BamTools::RefVector &getReferences() const { return _references; }

	const BamTools::BamAlignment &getAlignment() const { return _bamAlignment; }


	void getChrName(string &) const;
	int getBamChrId() const;
	int getStartPos() const;
	int getEndPos() const;
	void getName(string &) const;
	void getScore(string &) const;
	char getStrand() const;
	void getMateChrName(string &str) const;
	virtual int getNumFields() const { return MINIMUM_VALID_BAM_FIELDS; }

	refs_t* getCramRefs() { return _bamReader->GetReference(); }

protected:

	BamTools::BamReader *_bamReader;
	BamTools::BamAlignment _bamAlignment;
	bool _eof;
	string _bamHeader;
	BamTools::RefVector _references;
	bool _useTags;

	static const int MINIMUM_PRINTABLE_BAM_FIELDS = 6;
	static const int  MINIMUM_VALID_BAM_FIELDS = 11;
	void extractNameFromCore();
};


#endif /* BAMFILEREADER_H_ */
