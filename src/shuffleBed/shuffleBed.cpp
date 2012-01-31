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


BedShuffle::BedShuffle(string &bedFile, string &genomeFile, string &excludeFile, string &includeFile, 
                       bool haveSeed, bool haveExclude, bool haveInclude, bool sameChrom, 
                       float overlapFraction, int seed, bool chooseChrom) {

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


    // use the supplied seed for the random
    // number generation if given.  else,
    // roll our own.
    if (_haveSeed) {
        _seed = seed;
        srand(seed);
    }
    else {
        // thanks to Rob Long for the tip.
        _seed = (unsigned)time(0)+(unsigned)getpid();
        srand(_seed);
    }

    _bed         = new BedFile(bedFile);
    _genome      = new GenomeFile(genomeFile);
    _chroms      = _genome->getChromList();
    _numChroms   = _genome->getNumberOfChroms();
    _genomeSize  = _genome->getGenomeSize();

    if (_haveExclude) {
        _exclude = new BedFile(excludeFile);
        _exclude->loadBedFileIntoMap();
    }
    
    if (_haveInclude) {
        _include = new BedFile(includeFile);
        _include->loadBedFileIntoMapNoBin();
        
        _numIncludeChroms = 0;
        masterBedMapNoBin::const_iterator it    = _include->bedMapNoBin.begin(); 
        masterBedMapNoBin::const_iterator itEnd = _include->bedMapNoBin.end();
        for(; it != itEnd; ++it) {
            _includeChroms.push_back(it->first);
            _numIncludeChroms++;
        }
    }

    if (_haveExclude == true && _haveInclude == false)
        ShuffleWithExclusions();
    else if  (_haveExclude == false && _haveInclude == true)
        ShuffleWithInclusions();
    else
        Shuffle();
}


BedShuffle::~BedShuffle(void) {

}


void BedShuffle::Shuffle() {
    BED bedEntry;
    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            ChooseLocus(bedEntry);
            _bed->reportBedNewLine(bedEntry);
        }
    }
    _bed->Close();
}



void BedShuffle::ShuffleWithExclusions() {

    BED bedEntry;
    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            // keep looking as long as the chosen
            // locus happens to overlap with regions
            // that the user wishes to exclude.
            int  tries = 0;
            bool haveOverlap = false;
            do 
            {
                // choose a new locus
                ChooseLocus(bedEntry);
                haveOverlap = _exclude->anyHits(bedEntry.chrom, bedEntry.start, bedEntry.end,
                                                bedEntry.strand, false, false, _overlapFraction, false);
                tries++;
            } while ((haveOverlap == true) && (tries <= MAX_TRIES));
            

            if (tries > MAX_TRIES) {
                cerr << "Error, line " << _bed->_lineNum << ": tried " << MAX_TRIES << " potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on." << endl;
            }
            else {
                _bed->reportBedNewLine(bedEntry);
            }
        }
    }
    _bed->Close();
}


void BedShuffle::ShuffleWithInclusions() {

    BED bedEntry;     // used to store the current BED line from the BED file.

    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            // choose a new locus
            ChooseLocusFromInclusionFile(bedEntry);
            _bed->reportBedNewLine(bedEntry);
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
        do 
        {
            // we need to combine two consective calls to rand()
            // because RAND_MAX is 2^31 (2147483648), whereas
            // mammalian genomes are obviously much larger.
            uint32_t randStart = ((rand() << 31) | rand()) % _genomeSize;
            // use the above randomStart (e.g., for human 0..3.1billion) 
            // to identify the chrom and start on that chrom.
            pair<string, int> location = _genome->projectOnGenome(randStart);
            bedEntry.chrom = location.first;
            bedEntry.start = location.second;
            bedEntry.end   = bedEntry.start + length;
            chromSize      = _genome->getChromSize(location.first);
        } while (bedEntry.end > chromSize);
        // keep looking if we have exceeded the end of the chrom.
    }
    // OLD, quite arguably flawed, method.
    // 1. Choose a chrom randomly (i.e., not weighted by size)
    // 2. Choose a position on that chrom randomly
    else 
    {
        do 
        {
            if (_sameChrom == false) {
                randomChrom    = _chroms[rand() % _numChroms];
                chromSize      = _genome->getChromSize(randomChrom);
                randomStart    = rand() % chromSize;
                bedEntry.chrom = randomChrom;
                bedEntry.start = randomStart;
                bedEntry.end   = randomStart + length;
            }
            else {
                chromSize      = _genome->getChromSize(chrom);
                randomStart    = rand() % chromSize;
                bedEntry.start = randomStart;
                bedEntry.end   = randomStart + length;
            }
        } while (bedEntry.end > chromSize);
    }
}


void BedShuffle::ChooseLocusFromInclusionFile(BED &bedEntry) {

    string chrom    = bedEntry.chrom;
    CHRPOS length   = bedEntry.end - bedEntry.start;

    string randomChrom;
    CHRPOS randomStart;
    BED includeInterval;
    
    if (_sameChrom == false) {

        // grab a random chromosome from the inclusion file.
        randomChrom            = _includeChroms[rand() % _numIncludeChroms];
        // get the number of inclusion intervals for that chrom
        size_t size            =  _include->bedMapNoBin[randomChrom].size();
        // grab a random interval on the chosen chromosome.
        size_t interval        = rand() % size;
        // retreive a ranom -incl interval on the selected chrom
        includeInterval        = _include->bedMapNoBin[randomChrom][interval];

        bedEntry.chrom = randomChrom;
    }
    else {
        // get the number of inclusion intervals for the original chrom
        size_t size =  _include->bedMapNoBin[chrom].size();
        // grab a random interval on the chosen chromosome.
        includeInterval       = _include->bedMapNoBin[chrom][rand() % size];
    }
    
    randomStart    = includeInterval.start + rand() % (includeInterval.size());
    bedEntry.start = randomStart;
    bedEntry.end   = randomStart + length;
    
    // use recursion to ensure that the chosen location 
    // doesn't go past the end of the chrom
    if (bedEntry.end > ((size_t) _genome->getChromSize(chrom))) {
        //bedEntry.end = _genome->getChromSize(chrom);
        ChooseLocusFromInclusionFile(bedEntry);
    }
}

