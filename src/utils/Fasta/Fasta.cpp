// ***************************************************************************
// FastaIndex.cpp (c) 2010 Erik Garrison <erik.garrison@bc.edu>
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 9 February 2010 (EG)
// ---------------------------------------------------------------------------
// Modified to use htslib/faidx December 2022 (BP)

#include "Fasta.h"

using namespace std;

FastaIndexEntry::FastaIndexEntry(string name, CHRPOS length)
    : name(name)
    , length(length)
{}

FastaIndexEntry::FastaIndexEntry(void) // empty constructor
{ clear(); }

FastaIndexEntry::~FastaIndexEntry(void)
{}

void FastaIndexEntry::clear(void)
{
    name = "";
    length = 0;
}

FastaIndex::FastaIndex(bool useFullHeader) 
  : useFullHeader(useFullHeader)
{}

void FastaIndex::readIndexFile(string fname) {

    faidx = fai_load(fname.c_str());
    string idx_path = fname + ".fai";

    // NOTE: can use faidx_set_cache_size to improve performance of repeated queries.
    if(faidx == NULL) {
        cerr << "Warning: malformed fasta index file " << fname <<  endl;
        exit(1);
    }

    struct stat index_attr, fasta_attr;
    stat(idx_path.c_str(), &index_attr);
    stat(fname.c_str(), &fasta_attr);
    if(fasta_attr.st_mtime > index_attr.st_mtime) {
        cerr << "Warning: the index file is older than the FASTA file." << endl;
    }

    int n_seq = faidx_nseq(faidx);
    for(int i = 0; i < n_seq; i++){
        const char* chrom_name = faidx_iseq(faidx, i);
        int chrom_len = faidx_seq_len(faidx, chrom_name); // TODO: update to faidx_seq_len64 when htslib is updated.
        sequenceNames.push_back(chrom_name);
        this->insert(make_pair(chrom_name, FastaIndexEntry(chrom_name, chrom_len)));
    }
}

FastaIndex::~FastaIndex(void) {
    fai_destroy(faidx);
}

FastaIndexEntry FastaIndex::entry(string name) {
    FastaIndex::iterator e = this->find(name);
    if (e == this->end()) {
        cerr << "unable to find FASTA index entry for '" << name << "'" << endl;
        exit(1);
    } else {
        return e->second;
    }
}

bool FastaIndex::chromFound(string name) {
    FastaIndex::iterator e = this->find(name);
    if (e == this->end()) {
        return false;
    }
    return true;
}

string FastaIndex::indexFileExtension() { return ".fai"; }

void FastaReference::open(string reffilename, bool usemmap, bool useFullHeader) {
    filename = reffilename;

    if (useFullHeader)
      usingfullheader = true;
    index = new FastaIndex(useFullHeader);
    index->readIndexFile(reffilename);
}

FastaReference::~FastaReference(void) {
    delete index;
}

string FastaReference::getSequence(string seqname) {
    FastaIndexEntry entry = index->entry(seqname);

    int len = 0;
    // NOTE: could perhaps use string.reserve, then send the pointer to avoid the copy below.
    char *seq = fai_fetch(index->faidx, seqname.c_str(), &len); // TODO: update to fai_fetch64 when htslib is updated.
    if (len == -1) {
        cerr << "error getting sequence:" << seqname << " general error" << endl;
        exit(1);
    }
    if (len == -2) {
        cerr << "error getting sequence:" << seqname << " not found in index" << endl;
        exit(1);
    }

    string s = seq;
    std::free(seq);
    return s;
}

string FastaReference::getSubSequence(string seqname, CHRPOS start, CHRPOS length) {
    //cerr << "searching:" << seqname << ":" << start << "-" << start + length << endl;
    FastaIndexEntry entry = index->entry(seqname);
    if (start < 0 || length < 1) {
        cerr << "Error: cannot construct subsequence with negative offset or length < 1" << endl;
        exit(1);
    }
    int len = 0;
    char *seq = faidx_fetch_seq(index->faidx, seqname.c_str(), start, start + length - 1, &len); // TODO: update to fai_fetch_seq64 when htslib is updated.
    if (len == -1) {
        cerr << "error getting sequence:" << seqname << " general error" << endl;
        exit(1);
    }
    if (len == -2) {
        cerr << "error getting sequence:" << seqname << " not found in index" << endl;
        exit(1);
    }

    string s = seq;
    std::free(seq);
    return s;
}

CHRPOS FastaReference::sequenceLength(string seqname) {
    if (index->chromFound(seqname)) {
        FastaIndexEntry entry = index->entry(seqname);
        return entry.length;
    }
    // 0 length means the chrom wasn't found
    return 0;
}
