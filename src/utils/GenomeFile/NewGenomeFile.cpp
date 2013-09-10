/*****************************************************************************
  NewGenomeFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "NewGenomeFile.h"


NewGenomeFile::NewGenomeFile(const QuickString &genomeFilename)
: _maxId(-1)
{
    _genomeFileName = genomeFilename;
    loadGenomeFileIntoMap();
}

NewGenomeFile::NewGenomeFile(const BamTools::RefVector &refVector)
: _maxId(-1)
{
    for (size_t i = 0; i < refVector.size(); ++i) {
        QuickString chrom = refVector[i].RefName;
        CHRPOS length = refVector[i].RefLength;
        _maxId++;
        _chromSizeIds[chrom] = pair<CHRPOS, CHRPOS>(length, _maxId);
    }
}

// Destructor
NewGenomeFile::~NewGenomeFile(void) {
}


void NewGenomeFile::loadGenomeFileIntoMap() {

	FILE *fp = fopen(_genomeFileName.c_str(), "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: Can't open genome file %s. Exiting...\n", _genomeFileName.c_str());
		exit(1);
	}
	char sLine[2048];
	char chrName[2048];
	CHRPOS chrSize = 0;
	while (!feof(fp)) {
		memset(sLine, 0, 2048);
		memset(chrName, 0, 2048);
		chrSize = 0;
		fgets(sLine, 2048, fp);
		sscanf(sLine, "%s %d", chrName, &chrSize);
		if (strlen(sLine) == 0) {
			continue;
		}
		_maxId++;
		_chromSizeIds[chrName] = pair<CHRPOS, CHRPOS>(chrSize, _maxId);
		_startOffsets.push_back(_genomeLength);
		_genomeLength += chrSize;
		_chromList.push_back(chrName);

	}
	_startOffsets.push_back(_genomeLength); //insert the final length as the last element
	//to help with the lower_bound call in the projectOnGenome method.
	fclose(fp);
}

bool NewGenomeFile::projectOnGenome(CHRPOS genome_pos, QuickString &chrom, CHRPOS &start) {
    // search for the chrom that the position belongs on.
    // add 1 to genome position b/c of zero-based, half open.
    vector<CHRPOS>::const_iterator low =
        lower_bound(_startOffsets.begin(), _startOffsets.end(), genome_pos + 1);
    
    // use the iterator to identify the appropriate index 
    // into the chrom name and start vectors
    CHRPOS i = CHRPOS(low-_startOffsets.begin());
    if (i < 0 || i >= _chromList.size()) {
    	return false; //position not on genome
    }
    chrom = _chromList[i - 1];
    start = genome_pos - _startOffsets[i - 1];
    return true;
}
    
CHRPOS NewGenomeFile::getChromSize(const QuickString &chrom) {
	if (chrom == _currChromName) {
		return _currChromSize;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	_currChromName = iter->first;
    	_currChromSize = iter->second.first;
    	_currChromId = iter->second.second;
    	return _currChromSize;
    }
    return INT_MAX;
}

CHRPOS NewGenomeFile::getChromId(const QuickString &chrom) {
	if (chrom == _currChromName) {
		return _currChromId;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	_currChromName = iter->first;
    	_currChromSize = iter->second.first;
    	_currChromId = iter->second.second;
    	return _currChromId;
    }
    return INT_MAX;
}
