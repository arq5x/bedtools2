// ***************************************************************************
// BgzfStream_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011(DB)
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#include <api/internal/BgzfStream_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstring>
#include <algorithm>
using namespace std;

// constructor
BgzfStream::BgzfStream(void)
    : UncompressedBlockSize(Constants::BGZF_DEFAULT_BLOCK_SIZE)
    , CompressedBlockSize(Constants::BGZF_MAX_BLOCK_SIZE)
    , BlockLength(0)
    , BlockOffset(0)
    , BlockAddress(0)
    , IsOpen(false)
    , IsWriteOnly(false)
    , IsWriteCompressed(true)
    , Stream(NULL)
    , UncompressedBlock(NULL)
    , CompressedBlock(NULL)
{
    try {
        CompressedBlock   = new char[CompressedBlockSize];
        UncompressedBlock = new char[UncompressedBlockSize];
    } catch( std::bad_alloc& ba ) {
        fprintf(stderr, "BgzfStream ERROR: unable to allocate memory\n");
        exit(1);
    }
}

// destructor
BgzfStream::~BgzfStream(void) {
    if( CompressedBlock   ) delete[] CompressedBlock;
    if( UncompressedBlock ) delete[] UncompressedBlock;
}

// closes BGZF file
void BgzfStream::Close(void) {

    // skip if file not open
    if ( !IsOpen ) return;

    // if writing to file, flush the current BGZF block,
    // then write an empty block (as EOF marker)
    if ( IsWriteOnly ) {
        FlushBlock();
        int blockLength = DeflateBlock();
        fwrite(CompressedBlock, 1, blockLength, Stream);
    }

    // flush and close stream
    fflush(Stream);
    fclose(Stream);

    // reset flags
    IsWriteCompressed = true;
    IsOpen = false;
}

// compresses the current block
int BgzfStream::DeflateBlock(void) {

    // initialize the gzip header
    char* buffer = CompressedBlock;
    memset(buffer, 0, 18);
    buffer[0]  = Constants::GZIP_ID1;
    buffer[1]  = (char)Constants::GZIP_ID2;
    buffer[2]  = Constants::CM_DEFLATE;
    buffer[3]  = Constants::FLG_FEXTRA;
    buffer[9]  = (char)Constants::OS_UNKNOWN;
    buffer[10] = Constants::BGZF_XLEN;
    buffer[12] = Constants::BGZF_ID1;
    buffer[13] = Constants::BGZF_ID2;
    buffer[14] = Constants::BGZF_LEN;

    // set compression level
    const int compressionLevel = ( IsWriteCompressed ? Z_DEFAULT_COMPRESSION : 0 );

    // loop to retry for blocks that do not compress enough
    int inputLength = BlockOffset;
    int compressedLength = 0;
    unsigned int bufferSize = CompressedBlockSize;

    while ( true ) {

        // initialize zstream values
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)UncompressedBlock;
        zs.avail_in  = inputLength;
        zs.next_out  = (Bytef*)&buffer[Constants::BGZF_BLOCK_HEADER_LENGTH];
        zs.avail_out = bufferSize - Constants::BGZF_BLOCK_HEADER_LENGTH - Constants::BGZF_BLOCK_FOOTER_LENGTH;

        // initialize the zlib compression algorithm
        if ( deflateInit2(&zs,
                          compressionLevel,
                          Z_DEFLATED,
                          Constants::GZIP_WINDOW_BITS,
                          Constants::Z_DEFAULT_MEM_LEVEL,
                          Z_DEFAULT_STRATEGY) != Z_OK )
        {
            fprintf(stderr, "BgzfStream ERROR: zlib deflate initialization failed\n");
            exit(1);
        }

        // compress the data
        int status = deflate(&zs, Z_FINISH);
        if ( status != Z_STREAM_END ) {

            deflateEnd(&zs);

            // reduce the input length and try again
            if ( status == Z_OK ) {
                inputLength -= 1024;
                if ( inputLength < 0 ) {
                    fprintf(stderr, "BgzfStream ERROR: input reduction failed\n");
                    exit(1);
                }
                continue;
            }

            fprintf(stderr, "BgzfStream ERROR: zlib::deflateEnd() failed\n");
            exit(1);
        }

        // finalize the compression routine
        if ( deflateEnd(&zs) != Z_OK ) {
            fprintf(stderr, "BgzfStream ERROR: zlib::deflateEnd() failed\n");
            exit(1);
        }

        compressedLength = zs.total_out;
        compressedLength += Constants::BGZF_BLOCK_HEADER_LENGTH + Constants::BGZF_BLOCK_FOOTER_LENGTH;
        if ( compressedLength > Constants::BGZF_MAX_BLOCK_SIZE ) {
            fprintf(stderr, "BgzfStream ERROR: deflate overflow\n");
            exit(1);
        }

        break;
    }

    // store the compressed length
    BamTools::PackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

    // store the CRC32 checksum
    unsigned int crc = crc32(0, NULL, 0);
    crc = crc32(crc, (Bytef*)UncompressedBlock, inputLength);
    BamTools::PackUnsignedInt(&buffer[compressedLength - 8], crc);
    BamTools::PackUnsignedInt(&buffer[compressedLength - 4], inputLength);

    // ensure that we have less than a block of data left
    int remaining = BlockOffset - inputLength;
    if ( remaining > 0 ) {
        if ( remaining > inputLength ) {
            fprintf(stderr, "BgzfStream ERROR: after deflate, remainder too large\n");
            exit(1);
        }
        memcpy(UncompressedBlock, UncompressedBlock + inputLength, remaining);
    }

    // update block data
    BlockOffset = remaining;

    // return result
    return compressedLength;
}

// flushes the data in the BGZF block
void BgzfStream::FlushBlock(void) {

    // flush all of the remaining blocks
    while ( BlockOffset > 0 ) {

        // compress the data block
        int blockLength = DeflateBlock();

        // flush the data to our output stream
        int numBytesWritten = fwrite(CompressedBlock, 1, blockLength, Stream);
        if ( numBytesWritten != blockLength ) {
            fprintf(stderr, "BgzfStream ERROR: expected to write %u bytes during flushing, but wrote %u bytes\n",
                    blockLength, numBytesWritten);
            exit(1);
        }

        // update block data
        BlockAddress += blockLength;
    }
}

// decompresses the current block
int BgzfStream::InflateBlock(const int& blockLength) {

    // inflate the data from compressed buffer into uncompressed buffer
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)CompressedBlock + 18;
    zs.avail_in  = blockLength - 16;
    zs.next_out  = (Bytef*)UncompressedBlock;
    zs.avail_out = UncompressedBlockSize;

    int status = inflateInit2(&zs, Constants::GZIP_WINDOW_BITS);
    if ( status != Z_OK ) {
        fprintf(stderr, "BgzfStream ERROR: could not decompress block - zlib::inflateInit() failed\n");
        return -1;
    }

    status = inflate(&zs, Z_FINISH);
    if ( status != Z_STREAM_END ) {
        inflateEnd(&zs);
        fprintf(stderr, "BgzfStream ERROR: could not decompress block - zlib::inflate() failed\n");
        return -1;
    }

    status = inflateEnd(&zs);
    if ( status != Z_OK ) {
        fprintf(stderr, "BgzfStream ERROR: could not decompress block - zlib::inflateEnd() failed\n");
        return -1;
    }

    // return result
    return zs.total_out;
}

// opens the BGZF file for reading (mode is either "rb" for reading, or "wb" for writing)
bool BgzfStream::Open(const string& filename, const char* mode) {

    // close current stream, if necessary, before opening next
    if ( IsOpen ) Close();

    // determine open mode
    if ( strcmp(mode, "rb") == 0 )
        IsWriteOnly = false;
    else if ( strcmp(mode, "wb") == 0)
        IsWriteOnly = true;
    else {
        fprintf(stderr, "BgzfStream ERROR: unknown file mode: %s\n", mode);
        return false;
    }

    // open BGZF stream on a file
    if ( (filename != "stdin") && (filename != "stdout") )
        Stream = fopen(filename.c_str(), mode);

    // open BGZF stream on stdin
    else if ( (filename == "stdin") && (strcmp(mode, "rb") == 0 ) )
        Stream = freopen(NULL, mode, stdin);

    // open BGZF stream on stdout
    else if ( (filename == "stdout") && (strcmp(mode, "wb") == 0) )
        Stream = freopen(NULL, mode, stdout);

    if ( !Stream ) {
        fprintf(stderr, "BgzfStream ERROR: unable to open file %s\n", filename.c_str() );
        return false;
    }

    // set flag & return success
    IsOpen = true;
    return true;
}

// reads BGZF data into a byte buffer
int BgzfStream::Read(char* data, const unsigned int dataLength) {

    // if stream not open for reading (or empty request)
    if ( !IsOpen || IsWriteOnly || dataLength == 0 )
        return 0;

    // read blocks as needed until desired data length is retrieved
    char* output = data;
    unsigned int numBytesRead = 0;
    while ( numBytesRead < dataLength ) {

        // determine bytes available in current block
        int bytesAvailable = BlockLength - BlockOffset;

        // read (and decompress) next block if needed
        if ( bytesAvailable <= 0 ) {
            if ( !ReadBlock() ) return -1;
            bytesAvailable = BlockLength - BlockOffset;
            if ( bytesAvailable <= 0 ) break;
        }

        // copy data from uncompressed source buffer into data destination buffer
        char* buffer   = UncompressedBlock;
        int copyLength = min( (int)(dataLength-numBytesRead), bytesAvailable );
        memcpy(output, buffer + BlockOffset, copyLength);

        // update counters
        BlockOffset  += copyLength;
        output       += copyLength;
        numBytesRead += copyLength;
    }

    // update block data
    if ( BlockOffset == BlockLength ) {
        BlockAddress = ftell64(Stream);
        BlockOffset  = 0;
        BlockLength  = 0;
    }

    return numBytesRead;
}

// reads a BGZF block
bool BgzfStream::ReadBlock(void) {

    char header[Constants::BGZF_BLOCK_HEADER_LENGTH];
    int64_t blockAddress = ftell64(Stream);

    // read block header from file
    int count = fread(header, 1, sizeof(header), Stream);

    // if block header empty
    if ( count == 0 ) {
        BlockLength = 0;
        return true;
    }

    // if block header invalid size
    if ( count != sizeof(header) ) {
        fprintf(stderr, "BgzfStream ERROR: read block failed - could not read block header\n");
        return false;
    }

    // validate block header contents
    if ( !BgzfStream::CheckBlockHeader(header) ) {
        fprintf(stderr, "BgzfStream ERROR: read block failed - invalid block header\n");
        return false;
    }

    // copy header contents to compressed buffer
    int blockLength = BamTools::UnpackUnsignedShort(&header[16]) + 1;
    char* compressedBlock = CompressedBlock;
    memcpy(compressedBlock, header, Constants::BGZF_BLOCK_HEADER_LENGTH);
    int remaining = blockLength - Constants::BGZF_BLOCK_HEADER_LENGTH;

    // read remainder of block
    count = fread(&compressedBlock[Constants::BGZF_BLOCK_HEADER_LENGTH], 1, remaining, Stream);
    if ( count != remaining ) {
        fprintf(stderr, "BgzfStream ERROR: read block failed - could not read data from block\n");
        return false;
    }

    // decompress block data
    count = InflateBlock(blockLength);
    if ( count < 0 ) {
        fprintf(stderr, "BgzfStream ERROR: read block failed - could not decompress block data\n");
        return false;
    }

    // update block data
    if ( BlockLength != 0 )
        BlockOffset = 0;
    BlockAddress = blockAddress;
    BlockLength  = count;

    // return success
    return true;
}

// seek to position in BGZF file
bool BgzfStream::Seek(const int64_t& position) {

    // skip if not open
    if ( !IsOpen ) return false;

    // determine adjusted offset & address
    int     blockOffset  = (position & 0xFFFF);
    int64_t blockAddress = (position >> 16) & 0xFFFFFFFFFFFFLL;

    // attempt seek in file
    if ( fseek64(Stream, blockAddress, SEEK_SET) != 0 ) {
        fprintf(stderr, "BgzfStream ERROR: unable to seek in file\n");
        return false;
    }

    // update block data & return success
    BlockLength  = 0;
    BlockAddress = blockAddress;
    BlockOffset  = blockOffset;
    return true;
}

void BgzfStream::SetWriteCompressed(bool ok) {
    IsWriteCompressed = ok;
}

// get file position in BGZF file
int64_t BgzfStream::Tell(void) const {
    if ( !IsOpen )
        return 0;
    return ( (BlockAddress << 16) | (BlockOffset & 0xFFFF) );
}

// writes the supplied data into the BGZF buffer
unsigned int BgzfStream::Write(const char* data, const unsigned int dataLen) {

    // skip if file not open for writing
    if ( !IsOpen || !IsWriteOnly ) return false;

    // write blocks as needed til all data is written
    unsigned int numBytesWritten = 0;
    const char* input = data;
    unsigned int blockLength = UncompressedBlockSize;
    while ( numBytesWritten < dataLen ) {

        // copy data contents to uncompressed output buffer
        unsigned int copyLength = min(blockLength - BlockOffset, dataLen - numBytesWritten);
        char* buffer = UncompressedBlock;
        memcpy(buffer + BlockOffset, input, copyLength);

        // update counter
        BlockOffset     += copyLength;
        input           += copyLength;
        numBytesWritten += copyLength;

        // flush (& compress) output buffer when full
        if ( BlockOffset == blockLength ) FlushBlock();
    }

    // return result
    return numBytesWritten;
}
