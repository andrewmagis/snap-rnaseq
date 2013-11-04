/*++
Module Name:

    SAM.cpp

Abstract:

    Sequence Alignment Map (SAM) file writer and reader.

Environment:

    User mode service.

    SamWriter and SamReader (and their subclasses) aren't thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "SAM.h"
#include "Bam.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "Util.h"
#include "ReadSupplierQueue.h"
#include "FileFormat.h"
#include "AlignerOptions.h"
#include "directions.h"
#include "exit.h"
#include <vector>

using std::max;
using std::min;
using util::strnchr;

bool readIdsMatch(const char* id0, const char* id1)
{
    for (unsigned i = 0; ; i++) {
        char c0 = id0[i];
        char c1 = id1[i];

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == 0 || c0 == ' ' || c0 == '/') return true;  
    }
    return true;
}

bool readIdsMatch(Read *read0, Read *read1)
{
    if (read0->getIdLength() != read1->getIdLength()) {
        return false;
    }
    for (unsigned i = 0; i < read0->getIdLength(); i++) {
        char c0 = read0->getId()[i];
        char c1 = read1->getId()[i];

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == ' ' || c0 == '/') return true;  
    }
    return true;
}

    char *
strnchrs(char *str, char charToFind, char charToFind2, size_t maxLen) // Hokey version that looks for either of two chars
{
    for (size_t i = 0; i < maxLen; i++) {
        if (str[i] == charToFind || str[i] == charToFind2) {
            return str + i;
        }
        if (str[i] == 0) {
            return NULL;
        }
    }
    return NULL;
}

    char *
SAMReader::skipToBeyondNextRunOfSpacesAndTabs(char *str, const char *endOfBuffer, size_t *charsUntilFirstSpaceOrTab)
{
    if (NULL == str) return NULL;

    char *nextChar = str;
    while (nextChar < endOfBuffer && *nextChar != ' ' && *nextChar != '\n' && *nextChar != '\t' && *nextChar != '\r' /* for Windows CRLF text */) {
        nextChar++;
    }

    if (NULL != charsUntilFirstSpaceOrTab) {
        *charsUntilFirstSpaceOrTab = nextChar - str;
    }

    if (nextChar >= endOfBuffer || *nextChar == '\n') {
        return NULL;
    }

    while (nextChar < endOfBuffer && (' ' == *nextChar || '\t' == *nextChar || '\r' == *nextChar)) {
        nextChar++;
    }

    if (nextChar >= endOfBuffer) {
        return NULL;
    }

    return nextChar;
}


    SAMReader *
SAMReader::create(
    const DataSupplier* supplier,
    const char *fileName,
    const ReaderContext& context,
    _int64 startingOffset, 
    _int64 amountOfFileToProcess)
{
    DataReader* data = supplier->getDataReader(maxLineLen);
    SAMReader *reader = new SAMReader(data, context);
    reader->init(fileName, startingOffset, amountOfFileToProcess);
    return reader;
}

    void
SAMReader::readHeader(
    const char* fileName,
    ReaderContext& context)
{
    DataReader* data = DataSupplier::Default[false]->getDataReader(maxLineLen);
    if (! data->init(fileName)) {
        fprintf(stderr, "Unable to read file %s\n", fileName);
        soft_exit(1);
    }
    // todo: allow for larger headers
    _int64 headerSize = 1024 * 1024; // 1M header max
    char* buffer = data->readHeader(&headerSize);
    if (!parseHeader(fileName, buffer, buffer + headerSize, context.genome, &headerSize, &context.headerMatchesIndex)) {
        fprintf(stderr,"SAMReader: failed to parse header on '%s'\n",fileName);
        soft_exit(1);
    }
    _ASSERT(context.header == NULL);
    char* p = new char[headerSize + 1];
    memcpy(p, buffer, headerSize);
    p[headerSize] = 0;
    context.header = p;
    context.headerBytes = context.headerLength = headerSize;
    delete data;
}

SAMReader::SAMReader(
    DataReader* i_data,
    const ReaderContext& i_context)
    : ReadReader(i_context), data(i_data), headerSize(-1), clipping(i_context.clipping)
{
}


//
// Implement the ReadReader form of getNextRead, which doesn't include the
// alignment results by simply throwing them away.
//
    bool
SAMReader::getNextRead(Read *readToUpdate)
{
    return getNextRead(readToUpdate, NULL, NULL, NULL, NULL, NULL, NULL);
}

    bool
SAMReader::parseHeader(
    const char *fileName, 
    char *firstLine, 
    char *endOfBuffer, 
    const Genome *genome, 
    _int64 *o_headerSize,
    bool *o_headerMatchesIndex)
{
    char *nextLineToProcess = firstLine;
    *o_headerMatchesIndex = true;
    int numSQLines = 0;
    while (NULL != nextLineToProcess && nextLineToProcess < endOfBuffer && '@' == *nextLineToProcess) {
        if (!strncmp("@SQ",nextLineToProcess,3)) {
            //
            // These lines represent sequences in the reference genome, what are
            // called "pieces" in the Genome class.  (Roughly, chromosomes or major
            // variants like some versions of the MHC genes on chr6; or more
            // particularly the things that come in different FASTA files from the
            // reference assembly).
            //
            // Verify that they actually match what's in our reference genome.
            //
            numSQLines++;
            if (nextLineToProcess + 3 >= endOfBuffer || ' ' != nextLineToProcess[3] && '\t' != nextLineToProcess[3]) {
                fprintf(stderr,"Malformed SAM file '%s' has @SQ without a following space or tab.\n",fileName);
                return false;
            }

            char *snStart = nextLineToProcess + 4;
            while (snStart < endOfBuffer && strncmp(snStart,"SN:",__min(3,endOfBuffer-snStart)) && *snStart != '\n') {
                snStart++;
            }

            if (snStart >= endOfBuffer || *snStart == '\n') {
                fprintf(stderr,"Malformed @SQ line doesn't have 'SN:' in file '%s'\n",fileName);
                return false;
            }

            const size_t pieceNameBufferSize = 512;
            char pieceName[pieceNameBufferSize];
            for (unsigned i = 0; i < pieceNameBufferSize && snStart+3+i < endOfBuffer; i++) {
                if (snStart[3+i] == ' ' || snStart[3+i] == '\t' || snStart[3+i] == '\n') {
                    pieceName[i] = '\0';
                } else {
                    pieceName[i] = snStart[3+i];
                }
            }
            pieceName[pieceNameBufferSize - 1] = '\0';

            if (genome == NULL || !genome->getOffsetOfPiece(pieceName,NULL)) {
                *o_headerMatchesIndex = false;
            }
        } else if (!strncmp("@HD",nextLineToProcess,3) || !strncmp("@RG",nextLineToProcess,3) || !strncmp("@PG",nextLineToProcess,3) ||
            !strncmp("@CO",nextLineToProcess,3)) {
            //
            // Ignore these lines.
            //
        } else {
            fprintf(stderr,"Unrecognized header line in SAM file.\n");
            return false;
        }
        nextLineToProcess = strnchr(nextLineToProcess,'\n',endOfBuffer-nextLineToProcess) + 1;
    }

    *o_headerMatchesIndex &= genome != NULL && numSQLines == genome->getNumPieces();
    *o_headerSize = nextLineToProcess - firstLine;
    return true;
}

    bool
SAMReader::parseLine(char *line, char *endOfBuffer, char *result[], size_t *linelength, size_t fieldLengths[])
{
    *linelength = 0;

    char *next = line;
    char *endOfLine = strnchr(line,'\n',endOfBuffer-line);
    if (NULL == endOfLine) {
        return false;
    }

    //
    // Skip over any leading spaces and tabs
    //
    while (next < endOfLine && (*next == ' ' || *next == '\t')) {
        next++;
    }

    for (unsigned i = 0; i < nSAMFields; i++) {
        if (NULL == next || next >= endOfLine) {
            if (i == OPT) {
                // no optional fields
                result[OPT] = NULL;
                break;
            } else {
                //
                // Too few fields.
                //
                return false;
            }
        }

        result[i] = next;
        if (i == OPT) {
            // OPT field is actually all fields until end of line
            fieldLengths[OPT] = endOfLine - next;
            break;
        }

        next = skipToBeyondNextRunOfSpacesAndTabs(next,endOfLine,&fieldLengths[i]);
    }

    *linelength =  endOfLine - line + 1;    // +1 skips over the \n
    return true;
}

    void
SAMReader::getReadFromLine(
    const Genome        *genome,
    char                *line, 
    char                *endOfBuffer, 
    Read                *read, 
    AlignmentResult     *alignmentResult,
    unsigned            *genomeLocation, 
    Direction           *direction,
    unsigned            *mapQ,
    size_t              *lineLength,
    unsigned *           flag,
    const char **        cigar,
    ReadClippingType     clipping
    )
{
    char *field[nSAMFields];
    size_t fieldLength[nSAMFields];

    if (!parseLine(line, endOfBuffer, field, lineLength, fieldLength)) {
        fprintf(stderr, "Failed to parse SAM line:\n%.*s\n", lineLength, line);
        soft_exit(1);
    }

    //
    // We have to copy the piece name (RNAME) into its own buffer because the code in Genome expects
    // it to be a null-terminated string, while all we've got is one that's space delimited.
    //
    const size_t pieceNameBufferSize = 512;
    char pieceName[pieceNameBufferSize];
    unsigned offsetOfPiece;
    parsePieceName(genome, pieceName, pieceNameBufferSize, &offsetOfPiece, NULL, field, fieldLength);

    if (NULL != genomeLocation) {
        *genomeLocation = parseLocation(offsetOfPiece, field, fieldLength);
    }

    if (fieldLength[SEQ] != fieldLength[QUAL]) {
        fprintf(stderr,"SAMReader: QUAL string unequal in length to SEQ string.\n");
        soft_exit(1);
    }

    unsigned _flag;
    const size_t flagBufferSize = 20;   // More than enough
    char flagBuffer[flagBufferSize];
    if (fieldLength[FLAG] >= flagBufferSize) {
        fprintf(stderr,"SAMReader: flag field is too long.\n");
        soft_exit(1);
    }
    memcpy(flagBuffer,field[FLAG],fieldLength[FLAG]);
    flagBuffer[fieldLength[FLAG]] = '\0';
    if (1 != sscanf(flagBuffer,"%d",&_flag)) {
        fprintf(stderr,"SAMReader: couldn't parse FLAG field.\n");
        soft_exit(1);
    }

    if (NULL != read) {
        //
        // Clip reads where the quality strings end in '#'
        //
        read->init(field[QNAME],(unsigned)fieldLength[QNAME],field[SEQ],field[QUAL],(unsigned)fieldLength[SEQ]);
        //
        // If this read is RC in the SAM file, we need to reverse it here, since Reads are always the sense that they were as they came
        // out of the base caller.
        //

        if (_flag & SAM_REVERSE_COMPLEMENT) {
            read->becomeRC();
        }
        read->clip(clipping);

        if (field[OPT] != NULL) {
            unsigned n = (unsigned) fieldLength[OPT];
            while (n > 0 && (field[OPT][n-1] == '\n' || field[OPT][n-1] == '\r')) {
                n--;
            }
            read->setAuxiliaryData(field[OPT], n);
            for (char* p = field[OPT]; p != NULL && p < field[OPT] + fieldLength[OPT]; p = SAMReader::skipToBeyondNextRunOfSpacesAndTabs(p, field[OPT] + fieldLength[OPT])) {
                if (strncmp(p, "RG:Z:", 5) == 0) {
                    read->setReadGroup(READ_GROUP_FROM_AUX);
                    break;
                }
            }
        }
    }

    if (NULL != alignmentResult) {
        if (_flag & SAM_UNMAPPED) {
            *alignmentResult = NotFound;
        } else {
            if ('*' == pieceName[0]) {
                fprintf(stderr,"SAMReader: mapped read didn't have RNAME filled in.\n");
                soft_exit(1);
            }
            *alignmentResult = SingleHit;   // NB: This isn't quite right, we should look at MAPQ.
        }
    }

    if (NULL != direction) {
        *direction = (_flag & SAM_REVERSE_COMPLEMENT) ? RC : FORWARD;
    }

    if (NULL != mapQ) {
        *mapQ = atoi(field[MAPQ]);
        if (*mapQ > 255) {
            fprintf(stderr,"SAMReader: MAPQ field has bogus value\n");
            soft_exit(1);
        }
    }

    if (NULL != flag) {
        *flag = _flag;
    }

    if (NULL != cigar) {
        *cigar = field[CIGAR];
    }
}

    void
SAMReader::parsePieceName(
    const Genome* genome,
    char* pieceName,
    size_t pieceNameBufferSize,
    unsigned* o_offsetOfPiece,
	int* o_indexOfPiece,
    char* field[],
    size_t fieldLength[],
	unsigned rfield)
{
    if (fieldLength[rfield] >= pieceNameBufferSize) {  // >= because we need a byte for the \0
        fprintf(stderr,"SAMReader: too long an RNAME.  Can't parse.\n");
        soft_exit(1);
    }
    
    memcpy(pieceName,field[rfield],fieldLength[rfield]);
    pieceName[fieldLength[rfield]] = '\0';

    *o_offsetOfPiece = 0;
    if ('*' != pieceName[0] && genome != NULL && !genome->getOffsetOfPiece(pieceName,o_offsetOfPiece, o_indexOfPiece)) {
        //fprintf(stderr,"Unable to find piece '%s' in genome.  SAM file malformed.\n",pieceName);
        //soft_exit(1);
    }
}

    unsigned
SAMReader::parseLocation(
    unsigned offsetOfPiece,
    char* field[],
    size_t fieldLength[],
	unsigned rfield,
	unsigned posfield)
{
    unsigned oneBasedOffsetWithinPiece = 0;
    if ('*' != field[rfield][0] && '*' != field[posfield][0]) {
        //
        // We can't call sscanf directly into the mapped file, becuase it reads to the end of the
        // string even when it's satisfied all of its fields.  Since this can be gigabytes, it's not
        // really good for perf.  Instead, copy the POS field into a local buffer and null terminate it.
        //

        const unsigned posBufferSize = 20;
        char posBuffer[posBufferSize];
        if (fieldLength[posfield] >= posBufferSize) {
            fprintf(stderr,"SAMReader: POS field too long.\n");
            soft_exit(1);
        }
        memcpy(posBuffer,field[posfield],fieldLength[posfield]);
        posBuffer[fieldLength[posfield]] = '\0';
        if (0 == sscanf(posBuffer,"%d",&oneBasedOffsetWithinPiece)) {
            fprintf(stderr,"SAMReader: Unable to parse position when it was expected.\n");
            soft_exit(1);
        }
        if (0 == oneBasedOffsetWithinPiece) {
            fprintf(stderr,"SAMReader: Position parsed as 0 when it was expected.\n");
            soft_exit(1);
        }
        return offsetOfPiece + oneBasedOffsetWithinPiece - 1; // -1 is because our offset is 0 based, while SAM is 1 based.
    } else {
        return InvalidGenomeLocation;
    }
}
    
    void
SAMReader::init(
    const char *fileName,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    if (! data->init(fileName)) {
        fprintf(stderr, "Unable to read file %s\n", fileName);
        soft_exit(1);
    }
    headerSize = context.headerBytes;
    reinit(max(startingOffset, (_int64) context.headerBytes),
        amountOfFileToProcess == 0 || startingOffset >= (_int64) context.headerBytes ? amountOfFileToProcess
            : amountOfFileToProcess - (context.headerBytes - startingOffset));
}

    void
SAMReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
{
    _ASSERT(-1 != headerSize && startingOffset >= headerSize);  // Must call init() before reinit()
    data->reinit(startingOffset, amountOfFileToProcess);
    char* buffer;
    _int64 validBytes;
    if (!data->getData(&buffer, &validBytes)) {
        return;
    }
    if (startingOffset != headerSize) {
        char *firstNewline = strnchr(buffer,'\n',validBytes);
        if (NULL == firstNewline) {
            return;
        }

        data->advance((unsigned)(firstNewline - buffer + 1)); // +1 skips over the newline.
    }
}

    bool
SAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult,
    unsigned *genomeLocation,
    Direction *direction,
    unsigned *mapQ, 
    unsigned *flag,
    bool ignoreEndOfRange,
    const char **cigar)
{
    char* buffer;
    _int64 bytes;
    if (! data->getData(&buffer, &bytes)) {
        data->nextBatch();
        if (! data->getData(&buffer, &bytes)) {
            return false;
        }
    }
    char *newLine = strnchr(buffer, '\n', bytes);
    if (NULL == newLine) {
        //
        // There is no newline, so the line crosses the end of the buffer.
        // This should never happen since underlying reader manages overflow between chunks.
        //
        fprintf(stderr,"SAM file has too long a line, or doesn't end with a newline!  Failing.  fileOffset = %lld\n", data->getFileOffset());
        soft_exit(1);
    }

    size_t lineLength;
    read->setReadGroup(context.defaultReadGroup);
    getReadFromLine(context.genome, buffer,buffer + bytes, read, alignmentResult, genomeLocation, direction, mapQ, &lineLength, flag, cigar, clipping);
    read->setBatch(data->getBatch());
    data->advance((newLine + 1) - buffer);

    return true;
}

    ReadSupplierGenerator *
SAMReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context)
{
    //
    // single-ended SAM files always can be read with the range splitter.
    //
    RangeSplitter *splitter = new RangeSplitter(QueryFileSize(fileName), numThreads, 100);
    return new RangeSplittingReadSupplierGenerator(fileName, true, numThreads, context);
}
    
    PairedReadReader*
SAMReader::createPairedReader(
    const DataSupplier* supplier,
    const char *fileName,
    _int64 startingOffset,
    _int64 amountOfFileToProcess, 
    bool autoRelease,
    const ReaderContext& context)
{

    SAMReader* reader = SAMReader::create(DataSupplier::Default[false], fileName, context, 0, 0);
    if (reader == NULL) {
        return NULL;
    }
    return PairedReadReader::PairMatcher(reader, autoRelease);
}


    PairedReadSupplierGenerator *
SAMReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context)
{
    //
    // need to use a queue so that pairs can be matched
    //

    PairedReadReader* paired = SAMReader::createPairedReader(DataSupplier::Default[false], fileName, 0, 0, false, context);
    if (paired == NULL) {
        fprintf(stderr, "Cannot create reader on %s\n", fileName);
        soft_exit(1);
    }
    ReadSupplierQueue* queue = new ReadSupplierQueue(paired);
    queue->startReaders();
    return queue;
}

class SAMFormat : public FileFormat
{
public:
    SAMFormat(bool i_useM) : useM(i_useM) {}

    virtual bool isFormatOf(const char* filename) const;
    
    virtual void getSortInfo(const Genome* genome, char* buffer, _int64 bytes, unsigned* o_location, unsigned* o_readBytes, int* o_refID, int* o_pos) const;

    virtual ReadWriterSupplier* getWriterSupplier(AlignerOptions* options, const Genome* genome, const Genome* transcriptome, const GTFReader *gtf) const;

    virtual bool writeHeader(
        const ReaderContext& context, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine) const;

    virtual bool writeRead(
        const Genome * genome, const Genome * transcriptome, const GTFReader * gtf,
        LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace, 
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result, 
        int mapQuality, unsigned genomeLocation, Direction direction, bool isTranscriptome = false, unsigned tlocation = 0,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL, 
        AlignmentResult mateResult = NotFound, unsigned mateLocation = 0, Direction mateDirection = FORWARD,
        bool mateIsTranscriptome = false, unsigned mateTlocation = 0) const; 

private:
    static const char * computeCigarString(const Genome * genome, LandauVishkinWithCigar * lv,
        char * cigarBuf, int cigarBufLen, char * cigarBufWithClipping, int cigarBufWithClippingLen,
        const char * data, unsigned dataLength, unsigned basesClippedBefore, unsigned basesClippedAfter,
        unsigned genomeLocation, Direction direction, bool useM, int * editDistance, std::vector<unsigned> &tokens);

    const bool useM;
};

const FileFormat* FileFormat::SAM[] = { new SAMFormat(false), new SAMFormat(true) };

    bool
SAMFormat::isFormatOf(
    const char* filename) const
{
    return util::stringEndsWith(filename, ".sam");
}
    
    void
SAMFormat::getSortInfo(
    const Genome* genome,
    char* buffer,
    _int64 bytes,
	unsigned* o_location,
	unsigned* o_readBytes,
	int* o_refID,
	int* o_pos) const
{
    char* fields[SAMReader::nSAMFields];
    size_t lengths[SAMReader::nSAMFields];
    size_t lineLength;
    SAMReader::parseLine(buffer, buffer + bytes, fields, &lineLength, lengths);
    _ASSERT(lineLength < UINT32_MAX);
	if (o_readBytes != NULL) {
		*o_readBytes = (unsigned) lineLength;
	}
    if (lengths[SAMReader::POS] == 0 || fields[SAMReader::POS][0] == '*') {
		if (lengths[SAMReader::PNEXT] == 0 || fields[SAMReader::PNEXT][0] == '*') {
			if (o_location != NULL) {
				*o_location = UINT32_MAX;
			}
			if (o_refID != NULL) {
				*o_refID = -1;
			}
			if (o_pos != NULL) {
				*o_pos = 0;
			}
		} else {
			const size_t pieceNameBufferSize = 512;
			char pieceName[pieceNameBufferSize];
			unsigned offsetOfPiece;
			SAMReader::parsePieceName(genome, pieceName, pieceNameBufferSize, &offsetOfPiece, o_refID, fields, lengths, SAMReader::RNEXT);
			if (o_location != NULL) {
				*o_location = SAMReader::parseLocation(offsetOfPiece, fields, lengths, SAMReader::RNEXT, SAMReader::PNEXT);
			}
		}
    } else {
        const size_t pieceNameBufferSize = 512;
        char pieceName[pieceNameBufferSize];
        unsigned offsetOfPiece;
        SAMReader::parsePieceName(genome, pieceName, pieceNameBufferSize, &offsetOfPiece, o_refID, fields, lengths);
		if (o_location != NULL) {
	        *o_location = SAMReader::parseLocation(offsetOfPiece, fields, lengths);
		}
    }
}

    ReadWriterSupplier*
SAMFormat::getWriterSupplier(
    AlignerOptions* options,
    const Genome* genome, 
    const Genome* transcriptome,
    const GTFReader* gtf) const
{
    DataWriterSupplier* dataSupplier;
    if (options->sortOutput) {
        size_t len = strlen(options->outputFileTemplate);
        // todo: this is going to leak, but there's no easy way to free it, and it's small...
        char* tempFileName = (char*) malloc(5 + len);
        strcpy(tempFileName, options->outputFileTemplate);
        strcpy(tempFileName + len, ".tmp");
        dataSupplier = DataWriterSupplier::sorted(this, genome, tempFileName, options->sortMemory * (1ULL << 30),
            options->numThreads, options->outputFileTemplate, NULL);
    } else {
        dataSupplier = DataWriterSupplier::create(options->outputFileTemplate);
    }
    return ReadWriterSupplier::create(this, dataSupplier, genome, transcriptome, gtf);
}

    bool
SAMFormat::writeHeader(
    const ReaderContext& context,
    char *header,
    size_t headerBufferSize,
    size_t *headerActualSize,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine)
    const
{
    char *commandLine;
	size_t commandLineSize = 0;
	for (int i = 0; i < argc; i++) {
		commandLineSize += strlen(argv[i]) + 1;	// +1 is either a space or the terminating null
	}
	commandLine = new char[commandLineSize];
	commandLine[0] = '\0';
	for (int i = 0; i < argc; i++) {
		strcat(commandLine,argv[i]);
		if (i != argc-1) {
			strcat(commandLine," ");
		}
	}

    size_t bytesConsumed = snprintf(header, headerBufferSize, "@HD\tVN:1.4\tSO:%s\n%s%s@PG\tID:SNAP\tPN:SNAP\tCL:%s\tVN:%s\n", 
		sorted ? "coordinate" : "unsorted",
        context.header == NULL ? (rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine) : "",
        context.header == NULL ? "\n" : "",
        commandLine,version);

	delete [] commandLine;
	commandLine = NULL;
    if (bytesConsumed >= headerBufferSize) {
        fprintf(stderr,"SAMWriter: header buffer too small\n");
        return false;
    }

    if (context.header != NULL) {
		bool hasRG = false;
        for (const char* p = context.header; p < context.header + context.headerLength; ) {
            const char* newline = strnchr(p, '\n', (context.header + context.headerLength) - p);
            if (newline == NULL) {
                newline = context.header + context.headerLength;
            }
            _ASSERT(newline - p >= 3);
            // skip @HD lines, and also @SQ lines if header does not match index
            if (strncmp(p, "@HD", 3) != 0 && (context.headerMatchesIndex || strncmp(p, "@SQ", 3) != 0)) {
				hasRG |= strncmp(p, "@RG", 3) == 0;
                if (bytesConsumed + (newline - p) + 1 >= headerBufferSize) {
                    fprintf(stderr,"SAMWriter: header buffer too small\n");
                    return false;
                }
                memcpy(header + bytesConsumed, p, (newline - p));
                * (header + bytesConsumed + (newline - p)) = '\n';
                bytesConsumed += (newline - p) + 1;
            }
            p = newline + 1;
        }
		if (! hasRG) {
			int n = snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "%s\n",
				rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine);
			if (n > headerBufferSize - bytesConsumed) {
				fprintf(stderr, "SAMWriter: header buffer too small\n");
                return false;
            }
			bytesConsumed += n;
		}
    }
#ifndef SKIP_SQ_LINES
    if ((context.header == NULL || ! context.headerMatchesIndex) && context.genome != NULL) {
        // Write an @SQ line for each chromosome / piece in the genome
        const Genome::Piece *pieces = context.genome->getPieces();
        int numPieces = context.genome->getNumPieces();
        unsigned genomeLen = context.genome->getCountOfBases();
        for (int i = 0; i < numPieces; i++) {
            unsigned start = pieces[i].beginningOffset;
            unsigned end = (i + 1 < numPieces) ? pieces[i+1].beginningOffset : genomeLen;
            bytesConsumed += snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "@SQ\tSN:%s\tLN:%u\n", pieces[i].name, (end - start)-500);

            if (bytesConsumed >= headerBufferSize) {
                fprintf(stderr,"SAMWriter: header buffer too small\n");
                return false;
            }
        }
    }
#endif // SKIP_SQ_LINES

    *headerActualSize = bytesConsumed;
    return true;
}
    
    bool
getSAMData(
    const Genome * genome,
    const Genome * transcriptome,
    const GTFReader * gtf,
    LandauVishkinWithCigar * lv,
    // output data
    char* data,
    char* quality,
    unsigned dataSize,
    const char*& pieceName,
    int& pieceIndex,
    int& flags,
    unsigned& positionInPiece,
    int& mapQuality,
    const char*& matePieceName,
    int& matePieceIndex,
    unsigned& matePositionInPiece,
    _int64& templateLength,
    unsigned& fullLength,
    const char*& clippedData,
    unsigned& clippedLength,
    unsigned& basesClippedBefore,
    unsigned& basesClippedAfter,
    // input data
    size_t& qnameLen,
    Read * read,
    AlignmentResult result, 
    unsigned genomeLocation,
    Direction direction,
    bool isTranscriptome,
    bool useM,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    unsigned mateLocation,
    Direction mateDirection,
    bool mateIsTranscriptome)
{
    pieceName = "*";
    positionInPiece = 0;
    const char *cigar = "*";
    templateLength = 0;
    
    if (0 == qnameLen) {
         qnameLen = read->getIdLength();
    }
    
    //
    // If the aligner said it didn't find anything, treat it as such.  Sometimes it will emit the
    // best match that it found, even if it's not within the maximum edit distance limit (but will
    // then say NotFound).  Here, we force that to be SAM_UNMAPPED.
    //
    if (NotFound == result) {
        genomeLocation = InvalidGenomeLocation;
    }

    if (InvalidGenomeLocation == genomeLocation) {
        //
        // If it's unmapped, then always emit it in the forward direction.  This is necessary because we don't even include
        // the SAM_REVERSE_COMPLEMENT flag for unmapped reads, so there's no way to tell that we reversed it.
        //
        direction = FORWARD;
    }

    // Write the data and quality strings. If the read is reverse complemented, these need to
    // be backwards from the original read. Also, both need to be unclipped.
    clippedLength = read->getDataLength();
    fullLength = read->getUnclippedLength();
    if (fullLength > dataSize) {
        return false;
    }
    if (direction == RC) {
      for (unsigned i = 0; i < fullLength; i++) {
        data[fullLength - 1 - i] = COMPLEMENT[read->getUnclippedData()[i]];
        quality[fullLength - 1 - i] = read->getUnclippedQuality()[i];
      }
      clippedData = &data[fullLength - clippedLength - read->getFrontClippedLength()];
      basesClippedBefore = fullLength - clippedLength - read->getFrontClippedLength();
      basesClippedAfter = read->getFrontClippedLength();
    } else {
      memcpy(data, read->getUnclippedData(), read->getUnclippedLength());
      memcpy(quality, read->getUnclippedQuality(), read->getUnclippedLength());
      clippedData = read->getData();
      basesClippedBefore = read->getFrontClippedLength();
      basesClippedAfter = fullLength - clippedLength - basesClippedBefore;
    }

    int editDistance = -1;
    if (genomeLocation != InvalidGenomeLocation) {
        // This could be either a single hit read or a multiple hit read where we just
        // returned one location, but either way, let's print that location. We will then
        // set the quality to 60 if it was single or 0 if it was multiple. These are the
        // values the SAM FAQ suggests for aligners that don't compute confidence scores.
        if (direction == RC) {
            flags |= SAM_REVERSE_COMPLEMENT;
        }
        
        //genomeLocation has already been converted from transcriptome coordinates
        const Genome::Piece *piece = genome->getPieceAtLocation(genomeLocation);
        pieceName = piece->name;
        pieceIndex = (int)(piece - genome->getPieces());
        positionInPiece = genomeLocation - piece->beginningOffset + 1; // SAM is 1-based
        
        mapQuality = max(0, min(70, mapQuality));
    } else {
        flags |= SAM_UNMAPPED;
        mapQuality = 0;
    }

    if (hasMate) {
        flags |= SAM_MULTI_SEGMENT;
        flags |= (firstInPair ? SAM_FIRST_SEGMENT : SAM_LAST_SEGMENT);
        if (mateLocation != InvalidGenomeLocation) {
                    
            //genomeLocation has already been converted from transcriptome coordinates
            const Genome::Piece *matePiece = genome->getPieceAtLocation(mateLocation);
            matePieceName = matePiece->name;
            matePieceIndex = (int)(matePiece - genome->getPieces());
            matePositionInPiece = mateLocation - matePiece->beginningOffset + 1;

            if (mateDirection == RC) {
                flags |= SAM_NEXT_REVERSED;
            }

            if (genomeLocation == InvalidGenomeLocation) {
                //
                // The SAM spec says that for paired reads where exactly one end is unmapped that the unmapped
                // half should just have RNAME and POS copied from the mate.
                //
                pieceName = matePieceName;
                pieceIndex = matePieceIndex;
                matePieceName = "=";
                positionInPiece = matePositionInPiece;
            }

        } else {
            flags |= SAM_NEXT_UNMAPPED;
            //
            // The mate's unmapped, so point it at us.
            //
            matePieceName = "=";
            matePieceIndex = pieceIndex;
            matePositionInPiece = positionInPiece;
        }

        if (genomeLocation != InvalidGenomeLocation && mateLocation != InvalidGenomeLocation) {
            flags |= SAM_ALL_ALIGNED;
            // Also compute the length of the whole paired-end string whose ends we saw. This is slightly
            // tricky because (a) we may have clipped some bases before/after each end and (b) we need to
            // give a signed result based on whether our read is first or second in the pair.
            _int64 myStart = genomeLocation - basesClippedBefore;
            _int64 myEnd = genomeLocation + clippedLength + basesClippedAfter;
            _int64 mateBasesClippedBefore = mate->getFrontClippedLength();
            _int64 mateBasesClippedAfter = mate->getUnclippedLength() - mate->getDataLength() - mateBasesClippedBefore;
            _int64 mateStart = mateLocation - (mateDirection == RC ? mateBasesClippedAfter : mateBasesClippedBefore);
            _int64 mateEnd = mateLocation + mate->getDataLength() + (mateDirection == FORWARD ? mateBasesClippedAfter : mateBasesClippedBefore);
			if (pieceName == matePieceName) { // pointer (not value) comparison, but that's OK.
				if (myStart < mateStart) {
					templateLength = mateEnd - myStart;
				} else {
					templateLength = -(myEnd - mateStart);
				}
 			} // otherwise leave TLEN as zero.
        }

        if (pieceName == matePieceName) {
            matePieceName = "=";     // SAM Spec says to do this when they're equal (and not *, which won't happen because this is a pointer, not string, compare)
        }
    }
    return true;
}

    bool
SAMFormat::writeRead(
    const Genome * genome,
    const Genome * transcriptome,
    const GTFReader *gtf,
    LandauVishkinWithCigar * lv,
    char * buffer,
    size_t bufferSpace, 
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result, 
    int mapQuality,
    unsigned genomeLocation,
    Direction direction,
    bool isTranscriptome,
    unsigned tlocation,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    unsigned mateLocation,
    Direction mateDirection,
    bool mateIsTranscriptome,
    unsigned mateTlocation) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[cigarBufSize];
    char cigarNew[cigarBufSize];

    const int cigarBufWithClippingSize = MAX_READ * 2 + 32;
    char cigarBufWithClipping[cigarBufWithClippingSize];

    int flags = 0;
    const char *pieceName = "*";
    int pieceIndex = -1;
    unsigned positionInPiece = 0;
    const char *cigar = "*";
    const char *matePieceName = "*";
    int matePieceIndex = -1;
    unsigned matePositionInPiece = 0;
    _int64 templateLength = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore;
    unsigned basesClippedAfter;
    int editDistance = -1;
    
    if (! getSAMData(genome, transcriptome, gtf, lv, data, quality, MAX_READ, pieceName, pieceIndex, 
        flags, positionInPiece, mapQuality, matePieceName, matePieceIndex, matePositionInPiece, templateLength,
        fullLength, clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
        qnameLen, read, result, genomeLocation, direction, isTranscriptome, useM,
        hasMate, firstInPair, mate, mateResult, mateLocation, mateDirection, mateIsTranscriptome))
    {
        return false;
    }
    
    if (genomeLocation != InvalidGenomeLocation) {
    
        std::vector<unsigned> tokens;
        if (!isTranscriptome) {
    
            cigar = computeCigarString(genome, lv, cigarBuf, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize, 
                                       clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
                                       genomeLocation, direction, useM, &editDistance, tokens);
                                       
        } else {
        
            //Vector of tokens from CIGAR string for Transcriptome conversion
            cigar = computeCigarString(transcriptome, lv, cigarNew, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize, 
                                       clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
                                       tlocation, direction, useM, &editDistance, tokens);
                          
            //We need the pieceName for conversion             
            const Genome::Piece *transcriptomePiece = transcriptome->getPieceAtLocation(tlocation);
            const char* transcriptomePieceName = transcriptomePiece->name;
            unsigned transcriptomePositionInPiece = tlocation - transcriptomePiece->beginningOffset + 1; // SAM is 1-based
                                                                             
            //Insert splice junction
            //cigar = convertTranscriptomeToGenome(gtf, tokens, transcriptomePieceName, transcriptomePositionInPiece, cigarBufWithClipping, cigarBufWithClippingSize);
            lv->insertSpliceJunctions(gtf, tokens, transcriptomePieceName, transcriptomePositionInPiece, (char*) cigarBuf, cigarBufSize);
            cigar = cigarBuf;
            
        }
         
    }

    // Write the SAM entry, which requires the following fields:
    //
    // 1. QNAME: Query name of the read or the read pair
    // 2. FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
    // 3. RNAME: Reference sequence name
    // 4. POS: 1-Based leftmost position of clipped alignment
    // 5. MAPQ: Mapping quality (Phred-scaled)
    // 6. CIGAR: Extended CIGAR string (operations: MIDNSHP)
    // 7. MRNM: Mate reference name (‘=’ if same as RNAME)
    // 8. MPOS: 1-based leftmost mate position
    // 9. ISIZE: Inferred insert size
    // 10. SEQQuery: Sequence on the same strand as the reference
    // 11. QUAL: Query quality (ASCII-33=Phred base quality)    

    //
    // Some FASTQ files have spaces in their ID strings, which is illegal in SAM.  Just truncate them at the space.
    //
    const char *firstSpace = strnchr(read->getId(),' ',qnameLen);
    if (NULL != firstSpace) {
        qnameLen = (unsigned)(firstSpace - read->getId());
    }

    const int nmStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of editDistance
    char nmString[nmStringSize];  
    snprintf(nmString, nmStringSize, "\tNM:i:%d",editDistance);

    unsigned auxLen;
    bool auxSAM;
    char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
    static bool warningPrinted = false;
    const char* readGroupSeparator = "";
    const char* readGroupString = "";
    if (aux != NULL && (! auxSAM)) {
        if (! warningPrinted) {
            fprintf(stderr, "warning: translating optional fields from BAM->SAM not yet implemented, optional fields will not be included in output\n");
            warningPrinted = true;
        }
        if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
            for (BAMAlignAux* bamAux = (BAMAlignAux*) aux; (char*) bamAux < aux + auxLen; bamAux = bamAux->next()) {
                if (bamAux->tag[0] == 'R' && bamAux->tag[1] == 'G' && bamAux->val_type == 'Z') {
                    readGroupSeparator = "\tRG:Z:";
                    readGroupString = (char*) bamAux->value();
                    break;
                }
            }
        }
        aux = NULL;
        auxLen = 0;
    }
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        readGroupSeparator = "\tRG:Z:";
        readGroupString = read->getReadGroup();
    }
    int charsInString = snprintf(buffer, bufferSpace, "%.*s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%.*s\t%.*s%s%.*s%s%s\tPG:Z:SNAP%s\n",
        qnameLen, read->getId(),
        flags,
        pieceName,
        positionInPiece,
        mapQuality,
        cigar,
        matePieceName,
        matePositionInPiece,
        templateLength,
        fullLength, data,
        fullLength, quality,
        aux != NULL ? "\t" : "", auxLen, aux != NULL ? aux : "",
        readGroupSeparator, readGroupString,
        nmString);

    if (charsInString > bufferSpace) {
        //
        // Out of buffer space.
        //
        return false;
    } else if (charsInString == bufferSpace) {
      buffer[bufferSpace-1] = '\n'; // overwrite trailing null with newline
    }


    if (NULL != spaceUsed) {
        *spaceUsed = charsInString;
    }
    return true;
}

// Compute the CIGAR edit sequence string for a read against a given genome location.
// Returns this string if possible or "*" if we fail to compute it (which would likely
// be a bug due to lack of buffer space). The pointer returned may be to cigarBuf so it
// will only be valid until computeCigarString is called again.
    const char *
SAMFormat::computeCigarString(
    const Genome *              genome,
    LandauVishkinWithCigar *    lv,
    char *                      cigarBuf,
    int                         cigarBufLen,
    char *                      cigarBufWithClipping,
    int                         cigarBufWithClippingLen,
    const char *                data,
    unsigned                    dataLength,
    unsigned                    basesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    genomeLocation,
    Direction                   direction,
	bool						useM,
    int *                       editDistance,
    std::vector<unsigned>       &tokens
)
{
    const char *reference = genome->getSubstring(genomeLocation, dataLength);
    if (NULL != reference) {
        *editDistance = lv->computeEditDistance(
                            genome->getSubstring(genomeLocation, dataLength),
                            dataLength,
                            data,
                            dataLength,
                            MAX_K - 1,
                            cigarBuf,
                            cigarBufLen,
						    useM,
						    tokens);
    } else {
        //
        // Fell off the end of the chromosome.
        //
        return "*";
    }

    if (*editDistance == -2) {
        fprintf(stderr, "WARNING: computeEditDistance returned -2; cigarBuf may be too small\n");
        return "*";
    } else if (*editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            fprintf(stderr, "WARNING: computeEditDistance returned -1; this shouldn't happen\n");
            warningPrinted = true;
        }
        return "*";
    } else {
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        char clipBefore[16] = {'\0'};
        char clipAfter[16] = {'\0'};
        if (basesClippedBefore > 0) {
            snprintf(clipBefore, sizeof(clipBefore), "%uS", basesClippedBefore);
            tokens.insert(tokens.begin(), 'S');
            tokens.insert(tokens.begin(), basesClippedBefore);

        }
        if (basesClippedAfter > 0) {
            snprintf(clipAfter, sizeof(clipAfter), "%uS", basesClippedAfter);
            tokens.push_back(basesClippedAfter);
            tokens.push_back('S');
        }
        snprintf(cigarBufWithClipping, cigarBufWithClippingLen, "%s%s%s", clipBefore, cigarBuf, clipAfter);
        return cigarBufWithClipping;
    }
}

    
