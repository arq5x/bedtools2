#include "BamFileReader.h"
#include "ParseTools.h"
#include <cstdio>
BamFileReader::BamFileReader()
:  _bamReader(NULL),
   _eof(false),
   _useTags(true),
   _shouldDeleteBamReader(false)
{

}

BamFileReader::~BamFileReader()
{
	if (_bamReader != NULL && _shouldDeleteBamReader) {
		delete _bamReader;
		_shouldDeleteBamReader = false;
		_bamReader = NULL;
	}
}

bool BamFileReader::open()
{
//	if (_bamReader == NULL) {
//		_bamReader = new BamTools::BamReader();
//		_shouldDeleteBamReader = true;
//	}
//	if (_inputStream != NULL) {
//		try {
//			_bamReader->OpenStream(_inputStream);
//		}
//		catch (...) {
//			fprintf(stderr, "ERROR: Unable to open BAM file from standard input.\n");
//			exit(1);
//		}
//	} else {
//		try {
//			_bamReader->Open(_filename);
//		}
//		catch (...) {
//			fprintf(stderr, "ERROR: Unable to open BAM file %s\n", _filename.c_str());
//			exit(1);
//		}
//	}
    _bamHeader = _bamReader->GetHeaderText();
    _references = _bamReader->GetReferenceData();

	return true;
}

bool BamFileReader::isOpen() const
{
	return _bamReader->IsOpen();
}

void BamFileReader::close()
{
//	_bamReader->Close();
}

bool BamFileReader::readEntry()
{
	if (_useTags) {
		if (_bamReader->GetNextAlignment(_bamAlignment)) {
			return true;
		}
	} else {
		if (_bamReader->GetNextAlignmentCore(_bamAlignment)) {
			return true;
		}
	}
	//hit end of file
	_eof = true;
	return false;
}

void BamFileReader::getChrName(QuickString &str) const
{
	int refId = _bamAlignment.RefID;
	if (refId < 0) {
		return;
	}
	str =  _references[refId].RefName;
}

int BamFileReader::getBamChrId() const
{
	return _bamAlignment.RefID;
}

int BamFileReader::getStartPos() const
{
	return _bamAlignment.Position;
}

int BamFileReader::getEndPos() const
{
	return _bamAlignment.GetEndPosition(false, false);
}

void BamFileReader::getName(QuickString &str) const
{
	if (!_useTags) {
		str = _bamAlignment.SupportData.AllCharData.c_str();
	} else {
		str = _bamAlignment.Name;
	}
    if (_bamAlignment.IsFirstMate()) {
    	str += "/1";
    }
    else if (_bamAlignment.IsSecondMate()) {
    	str += "/2";
    }
}

void BamFileReader::getScore(QuickString &str) const
{
	int2str(_bamAlignment.MapQuality, str);
}

char BamFileReader::getStrand() const
{
	if (_bamAlignment.IsReverseStrand()) {
		return '-';
	}
	return '+';
}
