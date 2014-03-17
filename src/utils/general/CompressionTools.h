/*
 * CompressionTools.h
 *
 *  Created on: Apr 10, 2013
 *      Author: nek3d
 */

#ifndef COMPRESSIONTOOLS_H_
#define COMPRESSIONTOOLS_H_

#include <zlib.h>
#include "BTlist.h"

//pass an array of gzipped data, and an empty, PRE-ALLOCATED buffer to store the unzipped results.
//return value is the length of the uncompressed data.
unsigned long  inflateGzippedArray(const BTlist<int> &inbuf, //the input. Gzipped data.
		unsigned char* outbuf, //user must allocate enough space for uncompressed output!
		size_t outbufSize, //capacity of output buffer
		size_t inbufSize =0); //save time on call to strlen by providing size of inbuf.



#endif /* COMPRESSIONTOOLS_H_ */
