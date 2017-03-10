/*****************************************************************************
  shuffleBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "shuffleBed.h"


BedShuffle::BedShuffle(string &bedFile, string &genomeFile, 
                       string &excludeFile, string &includeFile, 
                       bool haveSeed, bool haveExclude,
                       bool haveInclude, bool sameChrom, 
                       float overlapFraction, int seed,
                       bool chooseChrom, bool isBedpe, size_t maxTries,
                       bool noOverlapping, bool preventExceedingChromEnd) 
{

    _bedFile         = bedFile;
    _genomeFile      = genomeFile;
    _excludeFile     = excludeFile;
    _includeFile     = includeFile;
    _sameChrom       = sameChrom;
    _haveExclude     = haveExclude;
    _haveInclude     = haveInclude;
    _overlapFraction = overlapFraction;
    _haveSeed        = haveSeed;
    _chooseChrom     = chooseChrom;
    _isBedpe         = isBedpe;
    _maxTries        = maxTries;
    _tries           = 0;
    _noOverlapping   = noOverlapping;
    _preventExceedingChromEnd = preventExceedingChromEnd;
    _cumLen          = 0;

    // use the supplied seed for the random
    // number generation if given.  else,
    // roll our own.
    if (_haveSeed) {
        _seed = seed;
    }
    else {
        // thanks to Rob Long for the tip.
        _seed = (unsigned)time(0)+(unsigned)getpid();
    }
#ifdef USE_RAND
    srand(_seed);
#else
    mt_rand.seed(_seed);
#endif
    
    if (_isBedpe == false)
        _bed         = new BedFile(bedFile);
    else
        _bedpe       = new BedFilePE(bedFile);

    _genome      = new GenomeFile(genomeFile);
    _chroms      = _genome->getChromList();
    _numChroms   = _genome->getNumberOfChroms();
    _genomeSize  = _genome->getGenomeSize();

    if (_haveExclude) {
        _exclude = new BedFile(excludeFile);
        _exclude->loadBedFileIntoMap();
    }
    else if (_noOverlapping) {
        // create an empty map that we add to as we iterate.
        _exclude = new BedFile();
        // force down correct code-path.
        _haveExclude = true;
    }
    
    if (_haveInclude) 
    {
        _include = new BedFile(includeFile);
        _include->loadBedFileIntoVector();
        _include->assignWeightsBasedOnSize();
        // for(std::vector<BED>::iterator it = _include->bedList.begin(); it != _include->bedList.end(); it++)
        // {
        //     _cumLen += (*it).end - (*it).start;
        // }	
    }

    if (_haveExclude == true && _haveInclude == false)
        ShuffleWithExclusions();
    else if  (_haveExclude == false && _haveInclude == true)
        ShuffleWithInclusions();
    else if (_haveExclude == true && _haveInclude == true)
        ShuffleWithInclusionsAndExclusions();
    else
        Shuffle();
}


BedShuffle::~BedShuffle(void) {

}


void BedShuffle::Shuffle() {
    
    if (_isBedpe == false) {
        BED bedEntry;
        _bed->Open();
        while (_bed->GetNextBed(bedEntry)) {
            if (_bed->_status == BED_VALID) {
                ChooseLocus(bedEntry);
                if(_noOverlapping){
                    _exclude->addBEDIntoMap(bedEntry);
                }
                if (_tries <= _maxTries)
                {
                    _bed->reportBedNewLine(bedEntry);                    
                }
            }
        }
        _bed->Close();
    }
    // BEDPE input
    else {
        int lineNum = 0;     // current input line number
        BedLineStatus status;

        BEDPE bedEntry, nullBedPE;
        _bedpe->Open();
        while ((status = 
            _bedpe->GetNextBedPE(bedEntry, lineNum)) != BED_INVALID) 
        {
            if (status == BED_VALID) {
                ChoosePairedLocus(bedEntry);
                _bedpe->reportBedPENewLine(bedEntry);
            }
            bedEntry = nullBedPE;
        }
        _bedpe->Close();
    }
}



void BedShuffle::ShuffleWithExclusions() {

    if (_isBedpe == false) {
        BED bedEntry;
        _bed->Open();
        while (_bed->GetNextBed(bedEntry)) {
            if (_bed->_status == BED_VALID) {
                // keep looking as long as the chosen
                // locus happens to overlap with regions
                // that the user wishes to exclude.
                _tries = 0;
                bool haveOverlap = false;
                do 
                {
                    // choose a new locus
                    ChooseLocus(bedEntry);
                    haveOverlap = _exclude->anyHits(bedEntry.chrom, 
                                                    bedEntry.start, 
                                                    bedEntry.end,
                                                    bedEntry.strand, 
                                                    false, false, 
                                                    _overlapFraction, false);
                    _tries++;
                } while ((haveOverlap == true) && (_tries <= _maxTries));


                if (_tries > _maxTries) {
                    cerr << "Error, line " << _bed->_lineNum 
                         << ": tried " << _maxTries 
                         << " potential loci for entry, but could not avoid "
                         << "excluded regions.  Ignoring entry and moving on." 
                         << endl;                }
                else {
                    if(_noOverlapping){
                        // future entries cannot overlap this one
                        _exclude->addBEDIntoMap(bedEntry);
                    }
                    _bed->reportBedNewLine(bedEntry);
                }
            }
        }
        _bed->Close();
    }
    // BEDPE input
    else {
        int lineNum = 0;     // current input line number
        BedLineStatus status;

        BEDPE bedEntry;
        _bedpe->Open();
        while ((status = 
                _bedpe->GetNextBedPE(bedEntry, lineNum)) != BED_INVALID) 
        {
            if (status == BED_VALID) {
                // keep looking as long as the chosen
                // locus happens to overlap with regions
                // that the user wishes to exclude.
                _tries = 0;
                bool haveOverlap1 = false;
                bool haveOverlap2 = false;
                do 
                {
                    // choose a new locus
                    ChoosePairedLocus(bedEntry);
                    haveOverlap1 = _exclude->anyHits(bedEntry.chrom1, 
                                                     bedEntry.start1, 
                                                     bedEntry.end1,
                                                     bedEntry.strand1, 
                                                     false, false, 
                                                     _overlapFraction, false);

                    haveOverlap2 = _exclude->anyHits(bedEntry.chrom2, 
                                                     bedEntry.start2, 
                                                     bedEntry.end2,
                                                     bedEntry.strand2, 
                                                     false, false, 
                                                     _overlapFraction, false);
                    _tries++;
                } while (((haveOverlap1 == true) || (haveOverlap2 == true))
                        && (_tries <= _maxTries));
                
                if (_tries > _maxTries) {
                    cerr << "Error, line " << _bed->_lineNum 
                         << ": tried " << _maxTries 
                         << " potential loci for entry, but could not avoid "
                         << "excluded regions.  Ignoring entry and moving on."
                         << endl;
                }
                else {
                    _bedpe->reportBedPENewLine(bedEntry);
                }
                
            }
        }
        _bedpe->Close();
    }
}


void BedShuffle::ShuffleWithInclusions() {

    BED bedEntry;     // used to store the current BED line from the BED file.
    CHRPOS chromSize;
    
    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            _tries = 0;
            // choose a new locus
            do {
                ChooseLocusFromInclusionFile(bedEntry);
                chromSize = _genome->getChromSize(bedEntry.chrom);
                _tries++;
            } while ((bedEntry.end > chromSize)
                    && (_tries <= _maxTries));
            if (_tries > _maxTries) {
                cerr << "Error, line " << _bed->_lineNum 
                     << ": tried " << _maxTries 
                     << " potential loci for entry, but could not avoid "
                     << "excluded regions.  Ignoring entry and moving on." 
                     << endl;                }
            else {
                _bed->reportBedNewLine(bedEntry);
           }
        }
    }
    _bed->Close();
}


void BedShuffle::ShuffleWithInclusionsAndExclusions() {

    BED bedEntry;     // used to store the current BED line from the BED file.

    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            
            // keep looking as long as the chosen
            // locus happens to overlap with regions
            // that the user wishes to exclude.
            _tries = 0;
            bool haveOverlap = false;
            do 
            {
                // choose a new locus
                ChooseLocusFromInclusionFile(bedEntry);
                haveOverlap = _exclude->anyHits(bedEntry.chrom, 
                                                bedEntry.start, 
                                                bedEntry.end,
                                                bedEntry.strand, 
                                                false, false, 
                                                _overlapFraction, false);
                _tries++;
            } while ((haveOverlap == true) && (_tries <= _maxTries));
            

            if (_tries > _maxTries) {
                cerr << "Error, line " << _bed->_lineNum 
                     << ": tried " << _maxTries 
                     << " potential loci for entry, but could not avoid "
                     << "excluded regions.  Ignoring entry and moving on." 
                     << endl;                }
            else {
                _bed->reportBedNewLine(bedEntry);
                if (_noOverlapping){
                    _exclude->addBEDIntoMap(bedEntry);
                }
            }
        }
    }
    _bed->Close();
}

void BedShuffle::ChooseLocus(BED &bedEntry) {

    string randomChrom;
    CHRPOS randomStart;
    CHRPOS chromSize;
    string chrom    = bedEntry.chrom;
    CHRPOS start    = bedEntry.start;
    CHRPOS end      = bedEntry.end;
    CHRPOS length   = end - start;

    // choose a position randomly among the _entire_ genome.
    if (_chooseChrom == false) 
    {
        _tries = 0;
        do 
        {

            // we need to combine two consecutive calls to rand()
            // because RAND_MAX is 2^31 (2147483648), whereas
            // mammalian genomes are obviously much larger.
#ifdef USE_RAND
            uint32_t randStart = ((((long) rand()) << 31) | rand()) %
                                 _genomeSize;
#else
            uint32_t randStart = ((((long) mt_rand()) << 31) | mt_rand()) %
                                 _genomeSize;
#endif
            // use the above randomStart (e.g., for human 0..3.1billion) 
            // to identify the chrom and start on that chrom.
            pair<string, int> location = _genome->projectOnGenome(randStart);
            bedEntry.chrom = location.first;
            bedEntry.start = location.second;
            bedEntry.end   = bedEntry.start + length;
            chromSize      = _genome->getChromSize(location.first);

            if ((bedEntry.end > chromSize) && 
                (_preventExceedingChromEnd == false))
            {
                bedEntry.end = chromSize;
                break;
            }
            _tries++;

        } while ((bedEntry.end > chromSize) && (_tries <= _maxTries));
        // keep looking if we have exceeded the end of the chrom.
        if (_tries > _maxTries) {
            cerr << "Error, line " << _bed->_lineNum 
                 << ": tried " << _maxTries 
                 << " potential loci for entry, but could not avoid "
                 << "excluded regions.  Ignoring entry and moving on." 
                 << endl;
        }
    }
    // OLD, quite arguably flawed, method.
    // 1. Choose a chrom randomly (i.e., not weighted by size)
    // 2. Choose a position on that chrom randomly
    else 
    {
        do 
        {
            if (_sameChrom == false) {
#ifdef USE_RAND
                randomChrom    = _chroms[rand() % _numChroms];
                chromSize      = _genome->getChromSize(randomChrom);
                randomStart    = rand() % chromSize;
#else
                randomChrom    = _chroms[mt_rand() % _numChroms];
                chromSize      = _genome->getChromSize(randomChrom);
                randomStart    = mt_rand() % chromSize;
#endif
                bedEntry.chrom = randomChrom;
                bedEntry.start = randomStart;
                bedEntry.end   = randomStart + length;
            }
            else {
                chromSize      = _genome->getChromSize(chrom);
#ifdef USE_RAND
                randomStart    = rand() % chromSize;
#else
                randomStart    = mt_rand() % chromSize;
#endif
                bedEntry.start = randomStart;
                bedEntry.end   = randomStart + length;
            }

            if ((bedEntry.end > chromSize) && 
                (_preventExceedingChromEnd == false))
            {
                bedEntry.end = chromSize;
                break;
            }

        } while (bedEntry.end > chromSize);
    }
}


void BedShuffle::ChoosePairedLocus(BEDPE &b) {
    
    CHRPOS foot1_len = b.end1 - b.start1;
    CHRPOS foot2_len = b.end2 - b.start2;
    CHRPOS length    = b.end2 - b.start1;

    if (b.chrom1 == b.chrom2)
    {
        CHRPOS chromSize;
        do 
        {
#ifdef USE_RAND
            uint32_t randStart = ((((long) rand()) << 31) | rand()) %
                                 _genomeSize;
#else
            uint32_t randStart = ((((long) mt_rand()) << 31) | mt_rand()) %
                                 _genomeSize;
#endif
            pair<string, int> location = _genome->projectOnGenome(randStart);
            b.chrom1  = location.first;
            b.chrom2  = location.first;
            b.start1  = location.second;
            b.end1    = b.start1 + foot1_len;
            b.end2    = b.start1 + length;
            b.start2  = b.end2 - foot2_len;
            chromSize      = _genome->getChromSize(location.first);
        } while ((b.end1 > chromSize) || 
                (b.start2 > chromSize) ||
                (b.end2 > chromSize));
        // keep looking if we have exceeded the end of the chrom.
    }
    else
    {
        CHRPOS chromSize1, chromSize2;
        do 
        {
#ifdef USE_RAND
            uint32_t rand1Start = ((((long) rand()) << 31) | rand()) %
                                  _genomeSize;
            uint32_t rand2Start = ((((long) rand()) << 31) | rand()) %
                                  _genomeSize;
#else
            uint32_t rand1Start = ((((long) mt_rand()) << 31) | mt_rand()) %
                                  _genomeSize;
            uint32_t rand2Start = ((((long) mt_rand()) << 31) | mt_rand()) %
                                  _genomeSize;
#endif
            pair<string, int> location1 = _genome->projectOnGenome(rand1Start);
            pair<string, int> location2 = _genome->projectOnGenome(rand2Start);
            
            b.chrom1  = location1.first;
            b.chrom2  = location2.first;
            b.start1  = location1.second;
            b.start2  = location2.second;
            
            b.end1    = b.start1 + foot1_len;
            b.end2    = b.start2 + foot2_len;
            chromSize1      = _genome->getChromSize(location1.first);
            chromSize2      = _genome->getChromSize(location2.first);
            
        } while ((b.end1 > chromSize1) || 
                (b.end2 > chromSize2));
        // keep looking if we have exceeded the end of the chrom.
    }
}


void BedShuffle::ChooseLocusFromInclusionFile(BED &bedEntry) 
{
    // choose an -incl interval randomly, yet weighted by the size of the incl interval.
    size_t length = (bedEntry.end - bedEntry.start);
#ifdef USE_RAND
    double runif =((double)rand()/(double) RAND_MAX);
#else
    double runif =((double)mt_rand()/(double) mt_rand.max());
#endif
    BED *includeInterval = _include->sizeWeightedSearch(runif);
    
    // choose a random start within the -incl interval and reconstruct shuffled record
#ifdef USE_RAND
    CHRPOS randomStart = includeInterval->start + (rand() % includeInterval->size());
#else
    CHRPOS randomStart = includeInterval->start + ((long) mt_rand() % includeInterval->size());
#endif
    bedEntry.chrom = includeInterval->chrom;
    bedEntry.start = randomStart;
    bedEntry.end   = randomStart + length;
}
