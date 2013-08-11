/*++

Module Name:

    bigalloc.h

Abstract:

    Headers for an allocator that uses big pages where appropriate and possible.

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

--*/

#pragma once

//#define PROFILE_BIGALLOC

#ifdef PROFILE_BIGALLOC
#define BigAlloc(s) BigAllocProfile((s), NULL, __FUNCTION__)

void *BigAllocProfile(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL,
        const char* caller = NULL);

#else
void *BigAlloc(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL);
#endif

void PrintBigAllocProfile();

void BigDealloc(void *memory);

void *BigReserve(
        size_t      sizeToReserve,
        size_t      *sizeReserved = NULL,
        size_t      *pageSize = NULL);

bool BigCommit(
        void        *memoryToCommit,
        size_t      sizeToCommit);

//
// This class is used to allocate a group of objects all onto a single set of big pages.  It requires knowing
// the amount of memory to be allocated when it's created.  It does not support deleting memory other than
// all at once.
//
class BigAllocator {
public:
    BigAllocator(size_t i_maxMemory);
    ~BigAllocator();

    virtual void *allocate(size_t amountToAllocate);
    virtual void assertAllMemoryUsed();

#if     _DEBUG
    void checkCanaries();
#else  // DEBUG
    void checkCanaries() {}
#endif  // DEBUG
private:

    char    *basePointer;
    char    *allocPointer;
    size_t  maxMemory;

#if     _DEBUG
    //
    // Stick a canary between each allocation and 
    unsigned    nCanaries;
    static const unsigned maxCanaries = 100;
    static const unsigned canaryValue = 0xca4a71e5;
    unsigned    *canaries[maxCanaries];
#endif  // DEBUG
};

//
// An allocator that doesn't actually allocate, it just counts how much it would allocate.  The idea is that
// you can write allocations in a fairly normal looking way, call them with this to see how much would be
// allocated, then create a real BigAllocator with that amount of memory.  That way, you don't need to
// keep in sync the actual allocation and the code that knows how much memory will be needed.
//
class CountingBigAllocator : public BigAllocator
{
public:
    CountingBigAllocator() :size(0), allocations(NULL), BigAllocator(0) {}
    ~CountingBigAllocator();

    virtual void *allocate(size_t amountToAllocate);
    virtual void assertAllMemoryUsed() {}
    size_t getMemoryUsed() {return size;}

private:
    size_t  size;

    struct Allocation {
        void *ptr;
        Allocation *next;
    } *allocations;
};

extern bool BigAllocUseHugePages;


// trivial per-thread heap for use in zalloc
struct ThreadHeap
{
    char* start;
    char* end;
    char* next;
    ThreadHeap(size_t bytes)
    {
        next = start = (char*) BigAlloc(bytes);
        end = start + bytes;
    }
    void* alloc(size_t bytes)
    {
        if (next + bytes <= end) {
            void* result = next;
            next += bytes;
            return result;
        }
        return NULL;
    }
    bool free(void* p)
    {
        return (char*)p >= start && (char*) p <= end;
    }
    void reset()
    {
        next = start;
    }
    ~ThreadHeap()
    {
        BigDealloc(start);
    }
};

void* zalloc(void* opaque, unsigned items, unsigned size);

void zfree(void* opaque, void* p);
