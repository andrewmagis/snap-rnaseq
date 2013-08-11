/*++

Module Name:

    BaseAligner.h

Abstract:

    Header for SNAP genome aligner

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "Aligner.h"
#include "LandauVishkin.h"
#include "BigAlloc.h"
#include "ProbabilityDistance.h"
#include "AlignerStats.h"
#include "directions.h"



class BaseAligner: public Aligner {
public:

    BaseAligner(
        GenomeIndex    *i_genomeIndex, 
        unsigned        i_maxHitsToConsider, 
        unsigned        i_maxK,
        unsigned        i_maxReadSize,
        unsigned        i_maxSeedsToUse,
        double          i_maxSeedCoverage,
        unsigned        i_extraSearchDepth,
        LandauVishkin<1>*i_landauVishkin = NULL,
        LandauVishkin<-1>*i_reverseLandauVishkin = NULL,
        AlignerStats   *i_stats = NULL,
        BigAllocator    *allocator = NULL);

    virtual ~BaseAligner();

        virtual AlignmentResult
    AlignRead(
        Read        *inputRead,
        unsigned    *genomeLocation,
        Direction   *hitDirection,
        int         *finalScore = NULL,
        int         *mapq = NULL);
        
    //
    // A richer version of AlignRead that allows for searching near a given location.
    // If searchRadius is not 0, constrain the search to distance maxSpacing around
    // searchLocation in the orientation given by searchRC.
    //
        AlignmentResult
    AlignRead(
        Read        *inputRead,
        unsigned    *genomeLocation,
        Direction   *hitDirection,
        int         *finalScore,
        int         *mapq,
        unsigned     searchRadius,       // If non-zero, constrain search around searchLocation in direction searchRC.
        unsigned     searchLocation,
        Direction    searchDirection);
        
    //
    // Statistics gathering.
    //

    _int64 getNHashTableLookups() const {return nHashTableLookups;}
    _int64 getLocationsScored() const {return nLocationsScored;}
    _int64 getNHitsIgnoredBecauseOfTooHighPopularity() const {return nHitsIgnoredBecauseOfTooHighPopularity;}
    _int64 getNReadsIgnoredBecauseOfTooManyNs() const {return nReadsIgnoredBecauseOfTooManyNs;}
    _int64 getNIndelsMerged() const {return nIndelsMerged;}
    void addIgnoredReads(_int64 newlyIgnoredReads) {nReadsIgnoredBecauseOfTooManyNs += newlyIgnoredReads;}

    const char *getRCTranslationTable() const {return rcTranslationTable;}

    inline int getMaxK() const {return maxK;}

    inline void setMaxK(int maxK_) {maxK = maxK_;}

    inline void setReadId(int readId_) {readId = readId_;}

    const char *getName() const {return "Base Aligner";}

    inline bool checkedAllSeeds() {return popularSeedsSkipped == 0;}

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(BaseAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/}
 
    inline bool getExplorePopularSeeds() {return explorePopularSeeds;}
    inline void setExplorePopularSeeds(bool newValue) {explorePopularSeeds = newValue;}

    inline bool getStopOnFirstHit() {return stopOnFirstHit;}
    inline void setStopOnFirstHit(bool newValue) {stopOnFirstHit = newValue;}

    static size_t getBigAllocatorReservation(bool ownLandauVishkin, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned numSeedsFromCommandLine, double seedCoverage);

private:

    bool hadBigAllocator;

    LandauVishkin<> *landauVishkin;
    LandauVishkin<-1> *reverseLandauVishkin;
    bool ownLandauVishkin;

    ProbabilityDistance *probDistance;

    // Maximum distance to merge candidates that differ in indels over.
    static const unsigned maxMergeDist = 48; // Must be even and <= 64

    char rcTranslationTable[256];

    _int64 nHashTableLookups;
    _int64 nLocationsScored;
    _int64 nHitsIgnoredBecauseOfTooHighPopularity;
    _int64 nReadsIgnoredBecauseOfTooManyNs;
    _int64 nIndelsMerged;

    //
    // A bitvector indexed by offset in the read indicating whether this seed is used.
    // This is here to avoid doing a memory allocation in the aligner.
    //
    BYTE *seedUsed;
    BYTE *seedUsedAsAllocated;  // Use this for deleting seedUsed.

    inline bool IsSeedUsed(unsigned indexInRead) const {
        return (seedUsed[indexInRead / 8] & (1 << (indexInRead % 8))) != 0;
    }

    inline void SetSeedUsed(unsigned indexInRead) {
        seedUsed[indexInRead / 8] |= (1 << (indexInRead % 8));
    }

    struct Candidate {
        Candidate() {init();}
        void init();

        unsigned        score;
        int             seedOffset;
    };

    static const unsigned hashTableElementSize = maxMergeDist;   // The code depends on this, don't change it

     struct HashTableElement {
        HashTableElement();
        void init();

        //
        // Doubly linked list for the weight buckets.
        //
        HashTableElement    *weightNext;
        HashTableElement    *weightPrev;

        //
        // Singly linked list for the hash table buckets.
        //
        HashTableElement    *next;

        _uint64             candidatesUsed;    // Really candidates we still need to score
        _uint64             candidatesScored;

        unsigned             baseGenomeLocation;
        unsigned             weight;
        unsigned             lowestPossibleScore;
        unsigned             bestScore;
        unsigned             bestScoreGenomeLocation;
        Direction            direction;
        bool                 allExtantCandidatesScored;
        double               matchProbabilityForBestScore;

        Candidate            candidates[hashTableElementSize];
    };

    //
    // Clearing out all of the pointers in the hash tables is expensive relative to running
    // an alignment, because usually the table is much bigger than the number of entries in it.
    // So, we avoid that expense by simply not clearing out the table at all.  Instead, along with
    // the pointers we keep an epoch number.  There's a corresponding epoch number in the
    // BaseAligner object, and if the two differ then the hash table bucket is empty.  We increment
    // the epoch number in the BaseAligner at the beginning of each alignment, thus effectively
    // clearing the hash table from the last run.
    //
    struct HashTableAnchor {
        HashTableElement *element;
        _int64            epoch;
    };

    _int64 hashTableEpoch;

    unsigned nUsedHashTableElements;
    unsigned hashTableElementPoolSize;
    HashTableElement *hashTableElementPool;

    const HashTableElement emptyHashTableElement;

    unsigned candidateHashTablesSize;
    HashTableAnchor *candidateHashTable[NUM_DIRECTIONS];
    
    HashTableElement *weightLists;
    unsigned highestUsedWeightList;

    static inline unsigned hash(unsigned key) {
        key = key * 131;    // Believe it or not, we spend a long time computing the hash, so we're better off with more table entries and a dopey function.
        return key;
    }

    static const unsigned UnusedScoreValue = 0xffff;

    // MAPQ parameters, currently not set to match Mason.  Using #define because VC won't allow "static const double".
#define SNP_PROB  0.001
#define GAP_OPEN_PROB  0.001
#define GAP_EXTEND_PROB  0.5

    //
    // Storage that's used during a call to AlignRead, but that's also needed by the
    // score function.  Since BaseAligner is single threaded, it's easier just to make
    // them member variables than to pass them around.
    //
    unsigned lowestPossibleScoreOfAnyUnseenLocation[NUM_DIRECTIONS];
    unsigned mostSeedsContainingAnyParticularBase[NUM_DIRECTIONS];
    unsigned nSeedsApplied[NUM_DIRECTIONS];
    unsigned bestScore;
    unsigned bestScoreGenomeLocation;
    unsigned secondBestScore;
    unsigned secondBestScoreGenomeLocation;
    int      secondBestScoreDirection;
    unsigned scoreLimit;
    unsigned lvScores;
    unsigned lvScoresAfterBestFound;
    double probabilityOfAllCandidates;
    double probabilityOfBestCandidate;
    int firstPassSeedsNotSkipped[NUM_DIRECTIONS];
    unsigned smallestSkippedSeed[NUM_DIRECTIONS];
    unsigned highestWeightListChecked;
    bool usedHammingThisAlignment;

    double totalProbabilityByDepth[AlignerStats::maxMaxHits];
    void updateProbabilityMass();

        bool
    score(
        bool             forceResult,
        Read            *read[NUM_DIRECTIONS],
        AlignmentResult *result,
        int             *finalScore,
        unsigned        *singleHitGenomeLocation,
        Direction       *hitDirection,
        int             *mapq);
    
    void clearCandidates();

    bool findElement(unsigned genomeLocation, Direction direction, HashTableElement **hashTableElement);
    void findCandidate(unsigned genomeLocation, Direction direction, Candidate **candidate, HashTableElement **hashTableElement);
    void allocateNewCandidate(unsigned genomeLoation, Direction direction, unsigned lowestPossibleScore, int seedOffset, Candidate **candidate, HashTableElement **hashTableElement);
    void incrementWeight(HashTableElement *element);
    void prefetchHashTableBucket(unsigned genomeLocation, Direction direction);

    const Genome *genome;
    GenomeIndex *genomeIndex;
    unsigned seedLen;
    unsigned maxHitsToConsider;
    unsigned maxK;
    unsigned maxReadSize;
    unsigned maxSeedsToUseFromCommandLine; // Max number of seeds to look up in the hash table
    double   maxSeedCoverage;  // Max seeds to used expressed as readSize/seedSize this is mutually exclusive with maxSeedsToUseFromCommandLine
    unsigned extraSearchDepth;
    unsigned numWeightLists;

    char *rcReadData;
    char *rcReadQuality;
    char *reversedRead[NUM_DIRECTIONS];

    unsigned nTable[256];

    int readId;
    
    // How many overly popular (> maxHits) seeds we skipped this run
    unsigned popularSeedsSkipped;

    bool explorePopularSeeds; // Whether we should explore the first maxHits hits even for overly
                              // popular seeds (useful for filtering reads that come from a database
                              // with many very similar sequences).

    bool stopOnFirstHit;      // Whether to stop the first time a location matches with less than
                              // maxK edit distance (useful when using SNAP for filtering only).

    AlignerStats *stats;
};
