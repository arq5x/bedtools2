/*
 * InflateStreamBuf.h
 *
 *  Created on: Apr 11, 2013
 *      Author: nek3d
 */

#ifndef INFLATESTREAMBUF_H_
#define INFLATESTREAMBUF_H_

using namespace std;

#include <iostream>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <zlib.h>
#include <cstdlib>

#define GZBUFSIZ BUFSIZ

class InflateStreamBuf:public std::streambuf
{
	void reset_decoder() {
		std::memset((void*)&strm,0,sizeof(z_stream));
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;
	}

public:
	InflateStreamBuf(std::istream* in):in(in),status_flag(0) {
		if(in==NULL) throw ("Null pointer");

		reset_decoder();

		bufsiz=GZBUFSIZ;
		buffer=(char*)std::malloc(GZBUFSIZ*sizeof(char));
		if(buffer==NULL) {
			throw ("out of memory");
		}
		if ( inflateInit2(&strm, 16+MAX_WBITS) != Z_OK) {
			throw ("::inflateInit failed");
		}
		setg(&buffer[0], &buffer[GZBUFSIZ], &buffer[GZBUFSIZ]);
	}

	virtual void close() {
		if(in!=NULL) {
			in=NULL;
			(void)inflateEnd(&strm);
			free(buffer);
			buffer=NULL;
		}
		in=NULL;
	}

	virtual ~InflateStreamBuf() {
		close();
	}

	virtual int underflow () {

		if(in==NULL) return EOF;
		if(status_flag == Z_STREAM_END) {
			close();
			return EOF;
		}
		in->read((char*)buffin, GZBUFSIZ);

		strm.avail_in = (uInt)in->gcount();

		if(strm.avail_in == 0) {
			close();
			return EOF;
		}
		strm.next_in=buffin;
	   /* run inflate() on input until output buffer not full */
		unsigned int total=0;
		do {
			strm.avail_out = GZBUFSIZ;
			strm.next_out = buffout;
			status_flag = ::inflate(&strm, Z_NO_FLUSH);

			switch (status_flag) {
				case 0:break;
				case Z_STREAM_END:
					   break;
				case Z_NEED_DICT:
					status_flag=Z_DATA_ERROR;
				case Z_DATA_ERROR:
				case Z_MEM_ERROR:
					close();
					throw status_flag;
					break;
				default:
					close();
					throw status_flag;
			}

			unsigned int have = GZBUFSIZ - strm.avail_out;
			if (total + have > bufsiz) {
				buffer=(char*)std::realloc(buffer,sizeof(char)*(total+have));
				if(buffer==NULL) throw ("out of memory");
				bufsiz = total + have;
			}
			memcpy((void*)&buffer[total], &buffout[0], sizeof(char)*have);
			total+=have;
			
			// Handle the concatenated GZIP files: It's part of gzip standard that a GZIP file can be concatenated directly, just like BGZF
			if(status_flag == Z_STREAM_END && (strm.avail_in != 0  || !in->eof())) {

				memmove(buffin, strm.next_in, strm.avail_in);
				int remaining = strm.avail_in;

				if(remaining > 0) {
					inflateEnd(&strm);
					reset_decoder();
					
					if ( inflateInit2(&strm, 16+MAX_WBITS) != Z_OK) {
						throw ("::inflateInit failed");
					}

					strm.avail_in = remaining;
					strm.next_in = buffin;
				}

				status_flag = 0;
			}
		} while (strm.avail_out == 0 && strm.avail_in > 0);



		setg(buffer, buffer, &buffer[total]);

		return total==0?EOF:this->buffer[0];
	}

protected:
	std::istream* in;
	unsigned char buffin[GZBUFSIZ];
	unsigned char buffout[GZBUFSIZ];
	char* buffer;
	size_t bufsiz;
	z_stream strm;
	int status_flag;
};



#endif /* INFLATESTREAMBUF_H_ */
