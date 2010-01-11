// ***************************************************************************
// BGZF.cpp (c) 2009 Derek Barnett, Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 11 January 2010 (DB)
// ---------------------------------------------------------------------------
// BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading & writing BGZF files
// ***************************************************************************

#include <algorithm>
#include "BGZF.h"
using namespace BamTools;
using std::string;
using std::min;

BgzfData::BgzfData(void)
    : UncompressedBlockSize(DEFAULT_BLOCK_SIZE)
    , CompressedBlockSize(MAX_BLOCK_SIZE)
    , BlockLength(0)
    , BlockOffset(0)
    , BlockAddress(0)
    , IsOpen(false)
    , IsWriteOnly(false)
    , Stream(NULL)
    , UncompressedBlock(NULL)
    , CompressedBlock(NULL)
{
    try {
        CompressedBlock   = new char[CompressedBlockSize];
        UncompressedBlock = new char[UncompressedBlockSize];
    } catch( std::bad_alloc& ba ) {
        printf("ERROR: Unable to allocate memory for our BGZF object.\n");
        exit(1);
    }
}

// destructor
BgzfData::~BgzfData(void) {
    if(CompressedBlock)   { delete[] CompressedBlock;   }
    if(UncompressedBlock) { delete[] UncompressedBlock; }
}

// closes BGZF file
void BgzfData::Close(void) {

	// skip if file not open, otherwise set flag
    if (!IsOpen) { return; }
    IsOpen = false;

    // flush the current BGZF block
    if (IsWriteOnly) { FlushBlock(); }

	// write an empty block (as EOF marker)
	int blockLength = DeflateBlock();
	fwrite(CompressedBlock, 1, blockLength, Stream);
	
    // flush and close
    fflush(Stream);
    fclose(Stream);
}

// compresses the current block
int BgzfData::DeflateBlock(void) {

    // initialize the gzip header
    char* buffer = CompressedBlock;
    memset(buffer, 0, 18);
    buffer[0]  = GZIP_ID1;
    buffer[1]  = (char)GZIP_ID2;
    buffer[2]  = CM_DEFLATE;
    buffer[3]  = FLG_FEXTRA;
    buffer[9]  = (char)OS_UNKNOWN;
    buffer[10] = BGZF_XLEN;
    buffer[12] = BGZF_ID1;
    buffer[13] = BGZF_ID2;
    buffer[14] = BGZF_LEN;

    // loop to retry for blocks that do not compress enough
    int inputLength = BlockOffset;
    int compressedLength = 0;
	unsigned int bufferSize = CompressedBlockSize;

    while(true) {
		
		// initialize zstream values
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)UncompressedBlock;
        zs.avail_in  = inputLength;
        zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
        zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
		
        // initialize the zlib compression algorithm
        if(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
            printf("ERROR: zlib deflate initialization failed.\n");
            exit(1);
        }

        // compress the data
        int status = deflate(&zs, Z_FINISH);
        if(status != Z_STREAM_END) {

            deflateEnd(&zs);

            // reduce the input length and try again
            if(status == Z_OK) {
                inputLength -= 1024;
                if(inputLength < 0) {
                    printf("ERROR: input reduction failed.\n");
                    exit(1);
                }
                continue;
            }

            printf("ERROR: zlib deflate failed.\n");
            exit(1);
        }

        // finalize the compression routine
        if(deflateEnd(&zs) != Z_OK) {
            printf("ERROR: deflate end failed.\n");
            exit(1);
        }

        compressedLength = zs.total_out;
        compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

        if(compressedLength > MAX_BLOCK_SIZE) {
            printf("ERROR: deflate overflow.\n");
            exit(1);
        }

        break;
    }

    // store the compressed length
    BgzfData::PackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

    // store the CRC32 checksum
    unsigned int crc = crc32(0, NULL, 0);
    crc = crc32(crc, (Bytef*)UncompressedBlock, inputLength);
    BgzfData::PackUnsignedInt(&buffer[compressedLength - 8], crc);
    BgzfData::PackUnsignedInt(&buffer[compressedLength - 4], inputLength);

    // ensure that we have less than a block of data left
    int remaining = BlockOffset - inputLength;
    if(remaining > 0) {
        if(remaining > inputLength) {
            printf("ERROR: remainder too large.\n");
            exit(1);
        }
        memcpy(UncompressedBlock, UncompressedBlock + inputLength, remaining);
    }

    BlockOffset = remaining;
    return compressedLength;
}

// flushes the data in the BGZF block
void BgzfData::FlushBlock(void) {

    // flush all of the remaining blocks
    while(BlockOffset > 0) {

        // compress the data block
        int blockLength = DeflateBlock();

        // flush the data to our output stream
        int numBytesWritten = fwrite(CompressedBlock, 1, blockLength, Stream);

        if(numBytesWritten != blockLength) {
            printf("ERROR: Expected to write %u bytes during flushing, but wrote %u bytes.\n", blockLength, numBytesWritten);
            exit(1);
		}
		
        BlockAddress += blockLength;
    }
}

// de-compresses the current block
int BgzfData::InflateBlock(const int& blockLength) {

    // Inflate the block in m_BGZF.CompressedBlock into m_BGZF.UncompressedBlock
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)CompressedBlock + 18;
    zs.avail_in  = blockLength - 16;
    zs.next_out  = (Bytef*)UncompressedBlock;
    zs.avail_out = UncompressedBlockSize;

    int status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK) {
        printf("inflateInit failed\n");
        exit(1);
    }

    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END) {
        inflateEnd(&zs);
        printf("inflate failed\n");
        exit(1);
    }

    status = inflateEnd(&zs);
    if (status != Z_OK) {
        printf("inflateEnd failed\n");
        exit(1);
    }

    return zs.total_out;
}

void BgzfData::Open(const string& filename, const char* mode) {

	// determine open mode
    if ( strcmp(mode, "rb") == 0 ) {
        IsWriteOnly = false;
    } else if ( strcmp(mode, "wb") == 0) {
        IsWriteOnly = true;
    } else {
        printf("ERROR: Unknown file mode: %s\n", mode);
        exit(1);
    }

	// open Stream to read to/write from file, stdin, or stdout
	// stdin/stdout option contributed by Aaron Quinlan (2010-Jan-03)
	if ( (filename != "stdin") && (filename != "stdout") ) {
		// read/wrtie BGZF data to/from a file
		Stream = fopen(filename.c_str(), mode);
	}
	else if ( (filename == "stdin") && (strcmp(mode, "rb") == 0 ) ) { 
		// read BGZF data from stdin
		Stream = freopen(NULL, mode, stdin);
	}
	else if ( (filename == "stdout") && (strcmp(mode, "wb") == 0) ) { 
		// write BGZF data to stdout
		Stream = freopen(NULL, mode, stdout);
	}
	
    if(!Stream) {
        printf("ERROR: Unable to open the BAM file %s\n", filename.c_str() );
        exit(1);
    }
    IsOpen = true;
}

int BgzfData::Read(char* data, const unsigned int dataLength) {

   if (dataLength == 0) { return 0; }

   char* output = data;
   unsigned int numBytesRead = 0;
   while (numBytesRead < dataLength) {

       int bytesAvailable = BlockLength - BlockOffset;
       if (bytesAvailable <= 0) {
           if ( ReadBlock() != 0 ) { return -1; }
           bytesAvailable = BlockLength - BlockOffset;
           if ( bytesAvailable <= 0 ) { break; }
       }

       char* buffer   = UncompressedBlock;
       int copyLength = min( (int)(dataLength-numBytesRead), bytesAvailable );
       memcpy(output, buffer + BlockOffset, copyLength);

       BlockOffset  += copyLength;
       output       += copyLength;
       numBytesRead += copyLength;
   }

   if ( BlockOffset == BlockLength ) {
       BlockAddress = ftell(Stream);
       BlockOffset  = 0;
       BlockLength  = 0;
   }

   return numBytesRead;
}

int BgzfData::ReadBlock(void) {

    char    header[BLOCK_HEADER_LENGTH];
    int64_t blockAddress = ftell(Stream);

    int count = fread(header, 1, sizeof(header), Stream);
    if (count == 0) {
        BlockLength = 0;
        return 0;
    }

    if (count != sizeof(header)) {
        printf("read block failed - count != sizeof(header)\n");
        return -1;
    }

    if (!BgzfData::CheckBlockHeader(header)) {
        printf("read block failed - CheckBlockHeader() returned false\n");
        return -1;
    }

    int blockLength = BgzfData::UnpackUnsignedShort(&header[16]) + 1;
    char* compressedBlock = CompressedBlock;
    memcpy(compressedBlock, header, BLOCK_HEADER_LENGTH);
    int remaining = blockLength - BLOCK_HEADER_LENGTH;

    count = fread(&compressedBlock[BLOCK_HEADER_LENGTH], 1, remaining, Stream);
    if (count != remaining) {
        printf("read block failed - count != remaining\n");
        return -1;
    }

    count = InflateBlock(blockLength);
    if (count < 0) { return -1; }

    if ( BlockLength != 0 ) {
        BlockOffset = 0;
    }

    BlockAddress = blockAddress;
    BlockLength  = count;
    return 0;
}

bool BgzfData::Seek(int64_t position) {

    int     blockOffset  = (position & 0xFFFF);
    int64_t blockAddress = (position >> 16) & 0xFFFFFFFFFFFFLL;

    if (fseek(Stream, blockAddress, SEEK_SET) != 0) {
        printf("ERROR: Unable to seek in BAM file\n");
        exit(1);
    }

    BlockLength  = 0;
    BlockAddress = blockAddress;
    BlockOffset  = blockOffset;
    return true;
}

int64_t BgzfData::Tell(void) {
    return ( (BlockAddress << 16) | (BlockOffset & 0xFFFF) );
}

// writes the supplied data into the BGZF buffer
unsigned int BgzfData::Write(const char* data, const unsigned int dataLen) {

    // initialize
    unsigned int numBytesWritten = 0;
    const char* input = data;
    unsigned int blockLength = UncompressedBlockSize;

    // copy the data to the buffer
    while(numBytesWritten < dataLen) {
        unsigned int copyLength = min(blockLength - BlockOffset, dataLen - numBytesWritten);
        char* buffer = UncompressedBlock;
        memcpy(buffer + BlockOffset, input, copyLength);

        BlockOffset     += copyLength;
        input           += copyLength;
        numBytesWritten += copyLength;

        if(BlockOffset == blockLength) {
            FlushBlock();
        }
    }

    return numBytesWritten;
}
