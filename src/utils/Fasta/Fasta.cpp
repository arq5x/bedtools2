// ***************************************************************************
// FastaIndex.cpp (c) 2010 Erik Garrison <erik.garrison@bc.edu>
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 9 February 2010 (EG)
// ---------------------------------------------------------------------------

#include "Fasta.h"

FastaIndexEntry::FastaIndexEntry(string name, int length, long long offset, int line_blen, int line_len)
    : name(name)
    , length(length)
    , offset(offset)
    , line_blen(line_blen)
    , line_len(line_len)
{}

FastaIndexEntry::FastaIndexEntry(void) // empty constructor
{ clear(); }

FastaIndexEntry::~FastaIndexEntry(void)
{}

void FastaIndexEntry::clear(void)
{
    name = "";
    length = NULL;
    offset = -1;  // no real offset will ever be below 0, so this allows us to
                  // check if we have already recorded a real offset
    line_blen = NULL;
    line_len  = NULL;
}

ostream& operator<<(ostream& output, const FastaIndexEntry& e) {
    // just write the first component of the name, for compliance with other tools
    output << split(e.name, ' ').at(0) << "\t" << e.length << "\t" << e.offset << "\t" <<
        e.line_blen << "\t" << e.line_len;
    return output;  // for multiple << operators.
}

FastaIndex::FastaIndex(void) 
{}

void FastaIndex::readIndexFile(string fname) {
    string line;
    long long linenum = 0;
    indexFile.open(fname.c_str(), ifstream::in);
    if (indexFile.is_open()) {
        while (getline (indexFile, line)) {
            ++linenum;
            // the fai format defined in samtools is tab-delimited, every line being:
            // fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen, (int)x.line_len
            vector<string> fields = split(line, '\t');
            if (fields.size() == 5) {  // if we don't get enough fields then there is a problem with the file
                // note that fields[0] is the sequence name
                char* end;
                string name = split(fields[0], " \t").at(0);  // key by first token of name
                sequenceNames.push_back(name);
                this->insert(make_pair(name, FastaIndexEntry(fields[0], atoi(fields[1].c_str()),
                                                    strtoll(fields[2].c_str(), &end, 10),
                                                    atoi(fields[3].c_str()),
                                                    atoi(fields[4].c_str()))));
            } else {
                cerr << "Warning: malformed fasta index file " << fname << 
                    "does not have enough fields @ line " << linenum << endl;
                cerr << line << endl;
                exit(1);
            }
        }
    } else {
        cerr << "could not open index file " << fname << endl;
        exit(1);
    }
}

// for consistency this should be a class method
bool fastaIndexEntryCompare ( FastaIndexEntry a, FastaIndexEntry b) { return (a.offset<b.offset); }

ostream& operator<<(ostream& output, FastaIndex& fastaIndex) {
    vector<FastaIndexEntry> sortedIndex;
    for(vector<string>::const_iterator it = fastaIndex.sequenceNames.begin(); it != fastaIndex.sequenceNames.end(); ++it)
    {
        sortedIndex.push_back(fastaIndex[*it]);
    }
    sort(sortedIndex.begin(), sortedIndex.end(), fastaIndexEntryCompare);
    for( vector<FastaIndexEntry>::iterator fit = sortedIndex.begin(); fit != sortedIndex.end(); ++fit) {
        output << *fit << endl;
    }
    return output;
}

void FastaIndex::indexReference(string refname) {
    // overview:
    //  for line in the reference fasta file
    //  track byte offset from the start of the file
    //  if line is a fasta header, take the name and dump the last sequnece to the index
    //  if line is a sequence, add it to the current sequence
    //cerr << "indexing fasta reference " << refname << endl;
    string line;
    FastaIndexEntry entry;  // an entry buffer used in processing
    entry.clear();
    int line_length = 0;
    long long offset = 0;  // byte offset from start of file
    long long line_number = 0; // current line number
    bool mismatchedLineLengths = false; // flag to indicate if our line length changes mid-file
                                        // this will be used to raise an error
                                        // if we have a line length change at
                                        // any line other than the last line in
                                        // the sequence
    bool emptyLine = false;  // flag to catch empty lines, which we allow for
                             // index generation only on the last line of the sequence
    ifstream refFile;
    refFile.open(refname.c_str());
    if (refFile.is_open()) {
        while (getline(refFile, line)) {
            ++line_number;
            line_length = line.length();
            if (line[0] == ';') {
                // fasta comment, skip
            } else if (line[0] == '+') {
                // fastq quality header
                getline(refFile, line);
                line_length = line.length();
                offset += line_length + 1;
                // get and don't handle the quality line
                getline(refFile, line);
                line_length = line.length();
            } else if (line[0] == '>' || line[0] == '@') { // fasta /fastq header
                // if we aren't on the first entry, push the last sequence into the index
                if (entry.name != "") {
                    mismatchedLineLengths = false; // reset line length error tracker for every new sequence
                    emptyLine = false;
                    flushEntryToIndex(entry);
                    entry.clear();
                }
                entry.name = line.substr(1, line_length - 1);
            } else { // we assume we have found a sequence line
                if (entry.offset == -1) // NB initially the offset is -1
                    entry.offset = offset;
                entry.length += line_length;
                if (entry.line_len) {
                    //entry.line_len = entry.line_len ? entry.line_len : line_length + 1;
                    if (mismatchedLineLengths || emptyLine) {
                        if (line_length == 0) {
                            emptyLine = true; // flag empty lines, raise error only if this is embedded in the sequence
                        } else {
                            if (emptyLine) {
                                cerr << "ERROR: embedded newline";
                            } else {
                                cerr << "ERROR: mismatched line lengths";
                            }
                            cerr << " at line " << line_number << " within sequence " << entry.name <<
                                endl << "File not suitable for fasta index generation." << endl;
                            exit(1);
                        }
                    }
                    // this flag is set here and checked on the next line
                    // because we may have reached the end of the sequence, in
                    // which case a mismatched line length is OK
                    if (entry.line_len != line_length + 1) {
                        mismatchedLineLengths = true;
                        if (line_length == 0) {
                            emptyLine = true; // flag empty lines, raise error only if this is embedded in the sequence
                        }
                    }
                } else {
                    entry.line_len = line_length + 1; // first line
                }
                entry.line_blen = entry.line_len - 1;
            }
            offset += line_length + 1;
        }
        // we've hit the end of the fasta file!
        // flush the last entry
        flushEntryToIndex(entry);
    } else {
        cerr << "could not open reference file " << refname << " for indexing!" << endl;
        exit(1);
    }
}

void FastaIndex::flushEntryToIndex(FastaIndexEntry& entry) {
    string name = split(entry.name, " \t").at(0);  // key by first token of name
    sequenceNames.push_back(name);
    this->insert(make_pair(name, FastaIndexEntry(entry.name, entry.length,
                        entry.offset, entry.line_blen,
                        entry.line_len)));

}

void FastaIndex::writeIndexFile(string fname) {
    //cerr << "writing fasta index file " << fname << endl;
    ofstream file;
    file.open(fname.c_str()); 
    if (file.is_open()) {
        file << *this;
    } else { 
        cerr << "could not open index file " << fname << " for writing!" << endl;
        exit(1);
    }
}

FastaIndex::~FastaIndex(void) {
    indexFile.close();
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

void FastaReference::open(string reffilename, bool usemmap) {
    filename = reffilename;
    if (!(file = fopen(filename.c_str(), "r"))) {
        cerr << "could not open " << filename << endl;
        exit(1);
    }
    index = new FastaIndex();
    struct stat stFileInfo; 
    string indexFileName = filename + index->indexFileExtension(); 
    // if we can find an index file, use it
    if(stat(indexFileName.c_str(), &stFileInfo) == 0) { 
        index->readIndexFile(indexFileName);
    } else { // otherwise, read the reference and generate the index file in the cwd
        cerr << "index file " << indexFileName << " not found, generating..." << endl;
        index->indexReference(filename);
        index->writeIndexFile(indexFileName);
    }
    if (usemmap) {
        usingmmap = true;
        int fd = fileno(file);
        struct stat sb;
        if (fstat(fd, &sb) == -1)
            cerr << "could not stat file" << filename << endl;
        filesize = sb.st_size;
        // map the whole file
        filemm = mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
    }
}

FastaReference::~FastaReference(void) {
    fclose(file);
    if (usingmmap) {
        munmap(filemm, filesize);
    }
    delete index;
}

string FastaReference::getSequence(string seqname) {
    FastaIndexEntry entry = index->entry(seqname);
    int newlines_in_sequence = entry.length / entry.line_blen;
    int seqlen = newlines_in_sequence  + entry.length;
    char* seq = (char*) calloc (seqlen + 1, sizeof(char));
    if (usingmmap) {
        memcpy(seq, (char*) filemm + entry.offset, seqlen);
    } else {
        fseek64(file, entry.offset, SEEK_SET);
        fread(seq, sizeof(char), seqlen, file);
    }
    seq[seqlen] = '\0';
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    pend = remove(pbegin, pend, '\n');
    pend = remove(pbegin, pend, '\0');
    string s = seq;
    free(seq);
    s.resize((pend - pbegin)/sizeof(char));
    return s;
}

// TODO cleanup; odd function.  use a map
string FastaReference::sequenceNameStartingWith(string seqnameStart) {
    try {
        return (*index)[seqnameStart].name;
    } catch (exception& e) {
        cerr << e.what() << ": unable to find index entry for " << seqnameStart << endl;
        exit(1);
    }
}

string FastaReference::getSubSequence(string seqname, int start, int length) {
    FastaIndexEntry entry = index->entry(seqname);
    if (start < 0 || length < 1) {
        cerr << "Error: cannot construct subsequence with negative offset or length < 1" << endl;
        exit(1);
    }
    // we have to handle newlines
    // approach: count newlines before start
    //           count newlines by end of read
    //             subtracting newlines before start find count of embedded newlines
    int newlines_before = start > 0 ? (start - 1) / entry.line_blen : 0;
    int newlines_by_end = (start + length - 1) / entry.line_blen;
    int newlines_inside = newlines_by_end - newlines_before;
    int seqlen = length + newlines_inside;
    char* seq = (char*) calloc (seqlen + 1, sizeof(char));
    if (usingmmap) {
        memcpy(seq, (char*) filemm + entry.offset + newlines_before + start, seqlen);
    } else {
        fseek64(file, (off_t) (entry.offset + newlines_before + start), SEEK_SET);
        fread(seq, sizeof(char), (off_t) seqlen, file);
    }
    seq[seqlen] = '\0';
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    pend = remove(pbegin, pend, '\n');
    pend = remove(pbegin, pend, '\0');
    string s = seq;
    free(seq);
    s.resize((pend - pbegin)/sizeof(char));
    return s;
}

long unsigned int FastaReference::sequenceLength(string seqname) {
    if (index->chromFound(seqname)) {
        FastaIndexEntry entry = index->entry(seqname);
        return entry.length;
    }
    // 0 length means the chrom wasn't found
    return 0;
}

