
#include "FileRecordTypeChecker.h"
#include "api/BamReader.h"
#include "ParseTools.h"

FileRecordTypeChecker::FileRecordTypeChecker()
: _eofHit(false),
  _inheader(false)
{
	_fileType = UNKNOWN_FILE_TYPE;
	_recordType = UNKNOWN_RECORD_TYPE;
	_numFields = 0;
	_isBinary = false;
	_isText = false;
	_isBed = false;
	_isDelimited = false;
	_delimChar = '\t'; //tab by default
	_firstValidDataLineIdx = -1;
	_isVCF = false;
	_isBAM = false;
	_isCRAM = false;
	_isGFF = false;
	_isGFFplus = false;
	_isGzipped = false;
	_isCompressed = false;
	_insufficientData = false;
	_fourthFieldNumeric = false;
	_givenEmptyBuffer = false;
	_isGroupBy = false;
	// TO DO: Bed4, Bed5, and BedGraph are missing from all of these.

	_hasName[UNKNOWN_RECORD_TYPE] = false;
	_hasName[EMPTY_RECORD_TYPE] = false;
	_hasName[BED3_RECORD_TYPE] = false;
	_hasName[BED6_RECORD_TYPE] = true;
	_hasName[BED12_RECORD_TYPE] = true;
	_hasName[BED_PLUS_RECORD_TYPE] = true;
	_hasName[BED6_PLUS_RECORD_TYPE] = true;
	_hasName[BAM_RECORD_TYPE] = true;
	_hasName[VCF_RECORD_TYPE] = true;
	_hasName[GFF_RECORD_TYPE] = true;
	_hasName[GFF_PLUS_RECORD_TYPE] = true;
	_hasName[NO_POS_PLUS_RECORD_TYPE] = true;

	_hasScore[UNKNOWN_RECORD_TYPE] = false;
	_hasScore[EMPTY_RECORD_TYPE] = false;
	_hasScore[BED3_RECORD_TYPE] = false;
	_hasScore[BED6_RECORD_TYPE] = true;
	_hasScore[BED12_RECORD_TYPE] = true;
	_hasScore[BED_PLUS_RECORD_TYPE] = true;
	_hasScore[BED6_PLUS_RECORD_TYPE] = true;
	_hasScore[BAM_RECORD_TYPE] = true;
	_hasScore[VCF_RECORD_TYPE] = true;
	_hasScore[GFF_RECORD_TYPE] = true;
	_hasScore[GFF_PLUS_RECORD_TYPE] = true;
	_hasScore[NO_POS_PLUS_RECORD_TYPE] = true;


	_hasStrand[UNKNOWN_RECORD_TYPE] = false;
	_hasStrand[EMPTY_RECORD_TYPE] = false;
	_hasStrand[BED3_RECORD_TYPE] = false;
	_hasStrand[BED6_RECORD_TYPE] = true;
	_hasStrand[BED12_RECORD_TYPE] = true;
	_hasStrand[BED_PLUS_RECORD_TYPE] = true; //actually, unknown. Give benefit of doubt.
	_hasStrand[BED6_PLUS_RECORD_TYPE] = true;
	_hasStrand[BAM_RECORD_TYPE] = true;
	_hasStrand[VCF_RECORD_TYPE] = true;
	_hasStrand[GFF_RECORD_TYPE] = true;
	_hasStrand[GFF_PLUS_RECORD_TYPE] = true;
	_hasStrand[NO_POS_PLUS_RECORD_TYPE] = true;


	_recordTypeNames[UNKNOWN_RECORD_TYPE] = "Unknown record type";
	_recordTypeNames[EMPTY_RECORD_TYPE] = "Empty record type";
	_recordTypeNames[BED3_RECORD_TYPE] = "Bed3 record type";
	_recordTypeNames[BED6_RECORD_TYPE] = "Bed6 record type";
	_recordTypeNames[BED12_RECORD_TYPE] = "Bed12 record type";
	_recordTypeNames[BED_PLUS_RECORD_TYPE] = "BedPlus record type";
	_recordTypeNames[BAM_RECORD_TYPE] = "BAM record type";
	_recordTypeNames[VCF_RECORD_TYPE] = "VCF record type";
	_recordTypeNames[GFF_RECORD_TYPE] = "Gff record type";
	_recordTypeNames[GFF_PLUS_RECORD_TYPE] = "GffPlus record type";
	_recordTypeNames[NO_POS_PLUS_RECORD_TYPE] = "NoPosPlus record type";


	_fileTypeNames[UNKNOWN_FILE_TYPE] = "Unknown file type";
	_fileTypeNames[EMPTY_FILE_TYPE] = "Empty file type";
	_fileTypeNames[SINGLE_LINE_DELIM_TEXT_FILE_TYPE] = "Delimited text file type";
	_fileTypeNames[GZIP_FILE_TYPE] = "Gzip file type";
	_fileTypeNames[BAM_FILE_TYPE] = "BAM file type";
	_fileTypeNames[VCF_FILE_TYPE] = "VCF file type";
}


bool FileRecordTypeChecker::scanBuffer(const char *buffer, size_t len, bool eofHit, bool isCompressed)
{
	_eofHit = eofHit;
	_isCompressed = isCompressed;
	_numBytesInBuffer = len;
	if (_numBytesInBuffer == 0) {
		_fileType = EMPTY_FILE_TYPE;
		_recordType = EMPTY_RECORD_TYPE;
		return true;
	}

	//special: the first thing we do is look for a gzipped file.
	if (!_isGzipped && ((unsigned char)(buffer[0]) == 0x1f)) {
		_isGzipped = true;
		return true;
	}
	//scan the first 8K block of the streamBuf.

	//now we have a buffer from the file.
	//first, test to see if it's binary or text.
	if (isBinaryBuffer(buffer, len)) {
		_isText = false;
		_isBinary = true;
		return true;
	} else {
		_isText = true;
		_isBinary = false;
		return handleTextFormat(buffer, len);
	}
}

bool FileRecordTypeChecker::isBinaryBuffer(const char *buffer, size_t len)
{
	if (isBAM(buffer)) {
		return true;
	}

	//Let's say that in a text file, at least 90% of the characters
	//should be alphanumeric, whitespace, or punctuation.
	static const float PERCENTAGE_PRINTABLE = .9;

	int alphaNumCount = 0;
	int whiteSpaceCount = 0;
	int punctuationCount = 0;

	for (int i=0; i < (int)len; i++) {
		char currChar = buffer[i];
		if (isalnum(currChar)) {
			alphaNumCount++;
		} else if (isspace(currChar)) {
			whiteSpaceCount++;
		} else if (ispunct(currChar)) {
			punctuationCount++;
		}
	}

	if ((float)(alphaNumCount + whiteSpaceCount + punctuationCount) / (float)(_numBytesInBuffer) < PERCENTAGE_PRINTABLE) {
		return true;
	}
	return false;
}


bool FileRecordTypeChecker::isBAM(const char *buffer)
{
	//check for BAM. The Bam Magic String is "BAM\1", and should be the first 4 characters of the file.

	if (strncmp(buffer, "BAM\1", 4) == 0) {
		_isBAM = true;
		_fileType = BAM_FILE_TYPE;
		_recordType = BAM_RECORD_TYPE;
		return true;
	}

	//TBD: Handle other binary formats
	return false;
}

bool FileRecordTypeChecker::handleTextFormat(const char *buffer, size_t len)
{
	if (isVCFformat(buffer)) {
		return isTextDelimtedFormat(buffer, len);
	} else if (isTextDelimtedFormat(buffer, len)) {
		//At this point, _isText and _isDelimited are set. _numFields and _delimChar are
		//set.
		_fileType = SINGLE_LINE_DELIM_TEXT_FILE_TYPE;
		if (_isGroupBy) {
			_recordType = NO_POS_PLUS_RECORD_TYPE;
			return true;
		}

		//Tokenize the first line of valid data into fields.
		//Need to make a copy so next call to tokenizer doesn't overwrite the line.

		string line(_tokenizer.getElem(_firstValidDataLineIdx));

		// ditch \r for Windows if necessary.
		if (line.size() && line[line.size()-1] == '\r') {
			line.resize(line.size()-1);
		}

		_tokenizer.setKeepFinalIncompleteElem(Tokenizer::USE_NOW);
		_tokenizer.setNumExpectedItems(_numFields);
		_tokenizer.tokenize(line, _delimChar);
		if (_tokenizer.getNumFields(line, _delimChar) != _numFields) {
			cerr << "Error: Type checker found wrong number of fields while tokenizing data line." << endl;
			cerr << "Perhaps you have extra TAB at the end of your line? Check with \"cat -t\""<< endl;
			exit(1);
		}

		if (isBedFormat()) {
			_isBed = true;
			if (_numFields == 3) {
				_recordType = BED3_RECORD_TYPE;
			} else if (_numFields == 4) {
				if (isNumeric(_tokenizer.getElem(3))) {
					_recordType = BEDGRAPH_RECORD_TYPE;
					_fourthFieldNumeric = true;
				} else {
					_fourthFieldNumeric = false;
					_recordType = BED4_RECORD_TYPE;
					_hasStrand[BED4_RECORD_TYPE] = isStrandField(3);
				}
			} else if (_numFields == 5 && passesBed5()) {
				_recordType = BED5_RECORD_TYPE;
			} else if (_numFields == 6 && passesBed6()) {
				_recordType = BED6_RECORD_TYPE;
			} else if (_numFields == 12 && passesBed12()) {
				_recordType = BED12_RECORD_TYPE;
			} else if (_numFields >3) {
				if (_numFields >= 6 && isStrandField(5)) {
					_recordType = BED6_PLUS_RECORD_TYPE;
				} else {
					_recordType = BED_PLUS_RECORD_TYPE;
				}

			}
			return true;
		}
		if (isGFFformat()) {
			if (_isGFFplus) {
				_recordType = GFF_PLUS_RECORD_TYPE;
				return true;
			}
			_isGFF = true;
			_recordType = GFF_RECORD_TYPE;
			return true;
		}
		//Here the Record must not have positions, so it is the NoPosPlus Type.
		_recordType = NO_POS_PLUS_RECORD_TYPE;
		return false;
	}
	return false;
}

bool FileRecordTypeChecker::isVCFformat(const char *buffer)
{
	if (_isVCF) {
		return true; //previous pass through this method has determined file is VCF.
	}
	if (memcmp(buffer, "##fileformat=VCF", 16) == 0) {
		_isVCF = true;
		_fileType = VCF_FILE_TYPE;
		_recordType = VCF_RECORD_TYPE;
		return true;
	}
	return false;
}

bool FileRecordTypeChecker::isBedFormat() {

	//test that the file has at least three fields.
	//2nd and 3rd fields of first valid data line must be integers. 3rd must not be less than 2nd.
	if (_numFields < 3) {
		return false;
	}
	//the 2nd and 3rd fields must be numeric.
	if (!isInteger(_tokenizer.getElem(1)) || !isInteger(_tokenizer.getElem(2))) {
		return false;
	}

	CHRPOS start = str2chrPos(_tokenizer.getElem(1));
	CHRPOS end = str2chrPos(_tokenizer.getElem(2));
	if (end < start) {
		return false;
	}
	return true;
}

bool FileRecordTypeChecker::isGFFformat()
{
	//a GFF file may have 8 or 9 fields. More than thats is GFFplus
	if (_numFields < 7 ) {
		return false;
	}
	//the 4th and 5th fields must be numeric.
	if (!isNumeric(_tokenizer.getElem(3)) || !isNumeric(_tokenizer.getElem(4))) {
		return false;
	}
	CHRPOS start = str2chrPos(_tokenizer.getElem(3));
	CHRPOS end = str2chrPos(_tokenizer.getElem(4));
	if (end < start) {
		return false;
	}
	if (_numFields > 8) {
		_isGFFplus = true;
	}
	return true;
}

bool FileRecordTypeChecker::isTextDelimtedFormat(const char *buffer, size_t len)
{
	//Break single string buffer into vector of strings. Delimiter is newline.
	_tokenizer.setKeepFinalIncompleteElem(Tokenizer::IGNORE);
	int numLines = _tokenizer.tokenize(buffer, '\n', _eofHit, _isCompressed);

	//anticipated delimiter characters are tab, comma, and semi-colon.
	//If we need new ones, they must be added in this method.
	//search each line for delimiter characters.

	vector<int> tabCounts;
	vector<int> commaCounts;
	vector<int> semicolonCounts;

	tabCounts.reserve(numLines);
	commaCounts.reserve(numLines);
	semicolonCounts.reserve(numLines);

	//loop through the lines, ignoring headers. Count potential delimiters,
	//see if we can find a few lines with the same number of a given delimiter.
	//delims are tested in hierarchical order, starting with tab,then comma, then semi-colon.

	int validLinesFound=0;
	int headerCount = 0;
	int emptyLines = 0;
	for (int i=0; i < numLines; i++ ) {


		if (validLinesFound >=4) {
			break; //really only need to look at like 4 lines of data, max.
		}

		const string line = _tokenizer.getElem(i);

		//skip over any empty line
		if (line.size() == 0) {
			emptyLines++;
			continue;
		}
		//
		//skip over any header line
		//

		if (_inheader) {
			headerCount++;
			_inheader = false; //inheaders can only apply to first line
			continue;
		}
		if (isHeaderLine(line)) {
			//clear any previously found supposedly valid data lines, because valid lines can only come after header lines.
			if (_firstValidDataLineIdx > -1 && _firstValidDataLineIdx < i) {
				_firstValidDataLineIdx = -1;
				validLinesFound--;
				headerCount++;
			}
			headerCount++;
			continue;
		}

		//a line must have some alphanumeric characters in order to be valid.
		bool hasAlphaNum = false;
		for (int j=0; j < len; j++) {
			if(isalnum(line[j])) {
				hasAlphaNum = true;
				break;
			}
		}
		if (!hasAlphaNum) {
			continue;
		}

		validLinesFound++;

		if (_firstValidDataLineIdx == -1) {
			_firstValidDataLineIdx = i;
		}

		int tab_count = std::count(line.begin(), line.end(), '\t');
		int comma_count = std::count(line.begin(), line.end(), ',');
		int semicolon_count = std::count(line.begin(), line.end(), ';');

		tabCounts.push_back(tab_count);
		commaCounts.push_back(comma_count);
		semicolonCounts.push_back(semicolon_count);
	}


	if (headerCount + emptyLines == numLines) {
		_insufficientData = true;
	}
	if (validLinesFound == 0) {
		return false;
	}
	_insufficientData = false;

	if (delimiterTesting(tabCounts, '\t')) {

		return true;
	}
	else if (validLinesFound) {
		return true;
	}

	if (delimiterTesting(commaCounts, ',')) {
		return true;
	}
	if (delimiterTesting(semicolonCounts, ';')) {
		return true;
	}

	return false; //unable to detect delimited file.
}

bool FileRecordTypeChecker::delimiterTesting(vector<int> &counts, char suspectChar)
{
	//check to see if we found the same number of tabs in every line.
	int numDelims = counts[0];
	if (numDelims != 0) {
		bool countsMatch = true;
		for (int i=1;  i < (int)counts.size(); i++) {
			if (counts[i] != numDelims) {
				countsMatch = false;
			}
		}
		if (countsMatch) {
			//Hurray!! We have successfully found a delimited file.
			_isDelimited = true;
			_delimChar = suspectChar;
			_numFields = numDelims + 1;
			return true;
		} else {
			return false;
		}
	}
	else { // there is just a single column with no delimiter.
		_numFields = 1;
		return false;
	}
}


void FileRecordTypeChecker::setBam()
{
	_fileType = BAM_FILE_TYPE;
	_recordType = BAM_RECORD_TYPE;
	_isBinary = true;
	_isBAM = true;
}

void FileRecordTypeChecker::setCram()
{
	_fileType = BAM_FILE_TYPE;
	_recordType = BAM_RECORD_TYPE;
	_isBinary = true;
	_isBAM = true;
	_isCRAM = true;
}

bool FileRecordTypeChecker::passesBed5() {
	return _isBed && _numFields == 5 && isNumeric(_tokenizer.getElem(4));
}

bool FileRecordTypeChecker::passesBed6() {
	return (_isBed && _numFields == 6 && isStrandField(5));
}

bool FileRecordTypeChecker::passesBed12() {

	return (isStrandField(5) && isNumeric(_tokenizer.getElem(6)) &&
			isNumeric(_tokenizer.getElem(7)) && isNumeric(_tokenizer.getElem(9)));
}

bool FileRecordTypeChecker::isStrandField(int field) {
	const string &strandChar = _tokenizer.getElem(field);
	return (strandChar == "+" || strandChar == "-" || strandChar == "." || strandChar == "*");
}
