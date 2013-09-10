/*
 * CompressionTools.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: nek3d
 */

#include "CompressionTools.h"
#include <cstring> //for strlen

unsigned long inflateGzippedArray(const BTlist<int> &inList, unsigned char* outbuf, size_t outbufSize, size_t inbufSize)
{
	if (inbufSize == 0) {
		inbufSize = inList.size();
	}
	Bytef *inbuf = new Bytef[inbufSize +1];
	memset(inbuf, 0, inbufSize +1);

	int i=0;
	for (BTlist<int>::const_iterator_type iter = inList.begin(); iter != inList.end(); iter = iter->next()) {
		inbuf[i] = iter->value();
		i++;
	}
	// zlib struct
	z_stream infstream;
	infstream.zalloc = Z_NULL;
	infstream.zfree = Z_NULL;
	infstream.opaque = Z_NULL;
	// setup "b" as the input and "c" as the compressed output
	infstream.avail_in = inbufSize; // size of input
	infstream.next_in = (Bytef *)inbuf; // input char array
	infstream.avail_out = outbufSize; // size of output
	infstream.next_out = (Bytef *)outbuf; // output char array

	// the actual DE-compression work.
	inflateInit2(&infstream, 16+MAX_WBITS); //decompresses gzipped data instead of zlib compressed data.
	inflate(&infstream, Z_NO_FLUSH);
	inflateEnd(&infstream);

	delete [] inbuf;

	return infstream.total_out;

}

