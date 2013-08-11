/*++

Module Name:

    DataWriter.h

Abstract:

    Headers for the DataWriter & related classes for the SNAP sequencer

Authors:

    Ravi Pandya, Feb 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "Compat.h"
#include "Read.h"

class DataWriterSupplier;

// per-thread writer for data into a single destination
class DataWriter
{
public:
    
    enum FilterType
    {
        ReadFilter, // reads data but does not modify it
        ModifyFilter, // modifies data in place
        CopyFilter, // copies data into new buffer, same size
        TransformFilter, // copies data into new buffer, possibly different size
    };
    // single filter instance per thread
    // points to filterSupplier for common data
    class Filter
    {
    public:
        Filter (FilterType i_filterType) : filterType(i_filterType) {}

        const FilterType filterType;

        virtual ~Filter() {}

		// called to set whether we're writing a header vs. individual reads
		virtual void inHeader(bool flag) {} // default do nothing

        // called when a chunk of data (i.e. a single read) has been written into the file
        virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, unsigned bytes, unsigned location) = 0;

        // called when a batch has been completed, after advancing to the next
        // e.g. so use getBatch(-1, ...) to get the one that was just completed
        // TransformFilters return #byte of transformed data in current buffer, so we need to advance again
        // TransformFilters should call getBatch(0) to ensure current buffer has been written before they write into it
        virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes) = 0;
    };
    
    // factory for per-thread filters
    class FilterSupplier
    {
    public:
        FilterSupplier (FilterType i_filterType) : filterType(i_filterType) {}
        
        const FilterType filterType;

        virtual ~FilterSupplier() {}

        FilterSupplier* compose(FilterSupplier* other);

        virtual Filter* getFilter() = 0;

        // called when entire file is done
        // TransformFilter will get a writer to append data if it wants
        virtual void onClose(DataWriterSupplier* supplier) = 0;
    };

    DataWriter(Filter* i_filter) : filter(i_filter) {}

    virtual ~DataWriter() {}

	void inHeader(bool flag)
	{ if (filter != NULL) { filter->inHeader(flag); } }

    // get remaining space in current buffer for writing
    virtual bool getBuffer(char** o_buffer, size_t* o_size) = 0;

    // advance within current buffer, reducing available space
    // should be called on each read, with the location
    virtual void advance(unsigned bytes, unsigned location = 0) = 0;

    // get complete data buffer in batch, relative==0 is current, relative==-1 is previous, etc.
    // if negative gets old data written, else waits for write to complete so you can write into it
    // o_offset gets physical offset (e.g. compressed), o_logical gets data offset (e.g. uncompressed)
    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size = NULL, size_t* o_used = NULL, size_t* o_offset = NULL, size_t* o_logicalUsed = 0, size_t* o_logicalOffset = NULL) = 0;

    // advance to next buffer
    virtual bool nextBatch() = 0;

    // this thread is complete
    virtual void close() = 0;

    // nanosecond timers
    static volatile _int64 FilterTime;
    static volatile _int64 WaitTime;

protected:
    Filter* filter;
};

class FileFormat;
class Genome;
class GzipWriterFilterSupplier;

// creates writers for multiple threads
class DataWriterSupplier
{
public:
    virtual ~DataWriterSupplier() {}

    virtual DataWriter* getWriter() = 0;

    // call when all threads are done, all filters destroyed
    virtual void close() = 0;
    
    static DataWriterSupplier* create(
        const char* filename,
        DataWriter::FilterSupplier* filterSupplier = NULL,
        int count = 4, size_t bufferSize = 16 * 1024 * 1024);
    
    static DataWriterSupplier* sorted(
        const FileFormat* format,
        const Genome* genome,
        const char* tempFileName,
        size_t tempBufferMemory,
        int numThreads,
        const char* sortedFileName,
        DataWriter::FilterSupplier* sortedFilterSupplier);

    // defaults follow BAM output spec
    static GzipWriterFilterSupplier* gzip(bool bamFormat, size_t chunkSize, int numThreads, bool bindToProcessors, bool multiThreaded);

    static DataWriter::FilterSupplier* markDuplicates(const Genome* genome);

    static DataWriter::FilterSupplier* bamIndex(const char* indexFileName, const Genome* genome, GzipWriterFilterSupplier* gzipSupplier);
};
