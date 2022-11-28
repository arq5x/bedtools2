// ***************************************************************************
// FastaIndex.h (c) 2010 Erik Garrison <erik.garrison@bc.edu>
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 February 2010 (EG)
// ---------------------------------------------------------------------------

#ifndef _FASTA_H
#define _FASTA_H

#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <algorithm>
#include "LargeFileSupport.h"
#include "bedFile.h"
#include "htslib/faidx.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include "split.h"
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

using namespace std;

class FastaIndexEntry {
    public:
        FastaIndexEntry(string name, CHRPOS length);
        FastaIndexEntry(void);
        ~FastaIndexEntry(void);
        string name;  // sequence name
        CHRPOS length;  // length of sequence
        void clear(void);
};

class FastaIndex : public map<string, FastaIndexEntry> {
    friend ostream& operator<<(ostream& output, FastaIndex& i);
    public:
        faidx_t *faidx;
        FastaIndex(bool useFullHeader);
        ~FastaIndex(void);
	    bool useFullHeader;
        vector<string> sequenceNames;
        void indexReference(string refName);
        void readIndexFile(string fname);
        FastaIndexEntry entry(string key);
        bool chromFound(string name);
        void flushEntryToIndex(FastaIndexEntry& entry);
        string indexFileExtension(void);
};

class FastaReference {
    public:
        void open(string reffilename, bool usemmap = false,
		bool useFullHeader = false);
        bool usingmmap;
        string filename;
        bool usingfullheader;
        FastaReference(void) : usingmmap(false), usingfullheader(false) { }
        ~FastaReference(void);
        size_t filesize;
        FastaIndex* index;
        vector<FastaIndexEntry> findSequencesStartingWith(string seqnameStart);
        string getSequence(string seqname);
        // potentially useful for performance, investigate
        // void getSequence(string seqname, string& sequence);
        string getSubSequence(string seqname, CHRPOS start, CHRPOS length);
        string sequenceNameStartingWith(string seqnameStart);
        CHRPOS sequenceLength(string seqname);
};

#endif
