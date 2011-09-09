/*****************************************************************************
  tagBam.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "tagBam.h"

// build
TagBam::TagBam(const string &bamFile, const vector<string> &annoFileNames,
            const vector<string> &annoLables, const string &tag,
            bool useNames, bool useScores, bool sameStrand, bool diffStrand, float overlapFraction):

    _bamFile(bamFile),
    _annoFileNames(annoFileNames),
    _annoLabels(annoLables),
    _tag(tag),
    _useNames(useNames),
    _useScores(useScores),
    _sameStrand(sameStrand),
    _diffStrand(diffStrand),
    _overlapFraction(overlapFraction)
{}


// destroy and delete the open file pointers
TagBam::~TagBam(void) {
    delete _bed;
    CloseAnnoFiles();
}


void TagBam::OpenAnnoFiles() {
    for (size_t i=0; i < _annoFileNames.size(); ++i) {
        BedFile *file = new BedFile(_annoFileNames[i]);
        file->loadBedFileIntoMap();
        _annoFiles.push_back(file);
    }
}


void TagBam::CloseAnnoFiles() {
    for (size_t i=0; i < _annoFiles.size(); ++i) {
        BedFile *file = _annoFiles[i];
        delete file;
        _annoFiles[i] = NULL;
    }
}


void TagBam::Tag() {

    // open the annotations files for processing;
    OpenAnnoFiles();

    // open the BAM file
    BamReader reader;
    BamWriter writer;
    reader.Open(_bamFile);
    // get header & reference information
    string bamHeader  = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // set compression mode
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
//    if ( _isUncompressedBam ) compressionMode = BamWriter::Uncompressed;
    writer.SetCompressionMode(compressionMode);
    // open our BAM writer
    writer.Open("stdout", bamHeader, refs);

    // rip through the BAM file and test for overlaps with each annotation file.
    BamAlignment al;
    vector<BED> hits;
    while (reader.GetNextAlignment(al)) {
        if (al.IsMapped() == true) {
            BED a;
            a.chrom = refs.at(al.RefID).RefName;
            a.start = al.Position;
            a.end   = al.GetEndPosition(false, false);
            a.strand = "+";
            if (al.IsReverseStrand()) a.strand = "-";
            
            ostringstream annotations;
            // annotate the BAM file based on overlaps with the annotation files.
            for (size_t i = 0; i < _annoFiles.size(); ++i) 
            {
                // grab the current annotation file.
                BedFile *anno = _annoFiles[i];
                
                if (!_useNames && !_useScores) {
                    // add the label for this annotation file to tag if there is overlap
                    if (anno->FindOneOrMoreOverlapsPerBin(a.chrom, a.start, a.end, a.strand, _sameStrand, _diffStrand, _overlapFraction)) 
                    {
                        annotations << _annoLabels[i] << ";";
                    }
                }
                // use the score field
                else if (!_useNames && _useScores) {
                    anno->FindOverlapsPerBin(a.chrom, a.start, a.end, a.strand, hits, _sameStrand, _diffStrand);
                    for (size_t i = 0; i < hits.size(); ++i) {
                        annotations << hits[i].score;
                        if (i < hits.size() - 1) annotations << ",";
                    }
                    if (hits.size() > 0) annotations << ";";
                    hits.clear();
                }
                // use the name field from the annotation files to populate tag
                else if (_useNames && !_useScores) {
                    anno->FindOverlapsPerBin(a.chrom, a.start, a.end, a.strand, hits, _sameStrand, _diffStrand);
                    for (size_t i = 0; i < hits.size(); ++i) {
                        annotations << hits[i].name;
                        if (i < hits.size() - 1) annotations << ",";
                    }
                    if (hits.size() > 0) annotations << ";";
                    hits.clear();
                }
            }
            // were there any overlaps with which to make a tag?
            if (annotations.str().size() > 0) {
                al.AddTag(_tag, "Z", annotations.str().substr(0, annotations.str().size() - 1)); // get rid of the last ";"
            }
            writer.SaveAlignment(al);
        }
    }
    reader.Close();

    // close the annotations files;
    CloseAnnoFiles();
}
