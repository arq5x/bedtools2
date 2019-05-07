/*
 * BamRecord.h
 *
 *  Created on: Dec 4, 2012
 *      Author: nek3d
 */

#ifndef BAMRECORD_H_
#define BAMRECORD_H_

#include "Bed6Interval.h"
#include "ParseTools.h"
#include "string.h"
#include "api/BamAlignment.h"

class FileReader;
class BamFileReader;
class RecordKeyVector;

class BamRecord : public Bed6Interval {
public:
	friend class FreeList<BamRecord>;

	BamRecord();
	virtual const BamRecord &operator=(const BamRecord &);


	// This using statement is only being added to supress warning from the CLANG compiler regarding
	// hidden overriden methods. Though it makes the base class methods available, developers should
	// not actually call them on a BamRecord object.
	using Bed6Interval::initFromFile;


	bool initFromFile(FileReader *);
	virtual bool initFromFile(BamFileReader *);
	virtual void clear();


	// As above, this using statement is only being added to supress warning from the CLANG compiler
	// regarding hidden overriden methods. Though it makes the base class methods available, developers
	// should not actually call them on a BamRecord object.
	using Bed6Interval::print;


	virtual void print(string &outBuf, int start, int end, RecordKeyVector *keyList) const;
	virtual void print(string &outBuf, RecordKeyVector *keyList) const;
	virtual void print(string &outBuf, const string & start, const string & end, RecordKeyVector *keyList) const;
	virtual void printNull(string &outBuf) const;
	virtual void printRemainingBamFields(string &outBuf, RecordKeyVector *keyList) const;
	virtual void printUnmapped(string &outBuf) const;

	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BAM_RECORD_TYPE; }
	const string &getCigarStr() const { return _cigarStr; }
	const vector<BamTools::CigarOp> &getCigarData() const { return _bamAlignment.CigarData; }

	const BamTools::BamAlignment &getAlignment() const { return _bamAlignment; }
	int getBamChromId() const { return _bamChromId; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 12; }
	static bool isNumericField(int fieldNum);

	CHRPOS getLength(bool obeySplits) const;

protected:
	BamTools::BamAlignment _bamAlignment;
	int _bamChromId; //different from chromId, because BAM file may be in different order
	//than the genomeFile.

	string _cigarStr; //stored for fast retrieval in column ops
	string _mateChrName;
	string _matePos;
	string _insertSize;
	string _queryBases;
	string _qualities;

	virtual ~BamRecord();
	void printRemainingBamFields();
	void buildCigarStr();

};


#endif /* BAMRECORD_H_ */
