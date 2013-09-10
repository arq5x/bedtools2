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
public:
	InflateStreamBuf(std::istream* in):in(in),status_flag(0) {
		if(in==NULL) throw ("Null pointer");
		std::memset((void*)&strm,0,sizeof(z_stream));
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;

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

		strm.avail_in = in->gcount();

		if(strm.avail_in == 0) {
			close();
			return EOF;
		}
		strm.next_in=buffin;
	   /* run inflate() on input until output buffer not full */
		unsigned int total=0;
		_currentOutBufSize = 0;
		do {
			strm.avail_out = GZBUFSIZ;
			strm.next_out = buffout;
			status_flag = ::inflate(&strm, Z_NO_FLUSH);

			switch (status_flag) {
				case 0:break;
				case Z_STREAM_END:break;
				case Z_NEED_DICT:
					status_flag=Z_DATA_ERROR;
				case Z_DATA_ERROR:
				case Z_MEM_ERROR:

					close();
					throw status_flag;
					break;
				default:
				{
				close();
				throw status_flag;
				break;
				}
			}

			unsigned int have = GZBUFSIZ - strm.avail_out;
			buffer=(char*)std::realloc(buffer,sizeof(char)*(total+have));
			if(buffer==NULL) {
				throw ("out of memory");
			}

			memcpy((void*)&buffer[total], &buffout[0], sizeof(char)*have);
			total+=have;
		} while (strm.avail_out == 0 && strm.avail_in > 0);


		setg(buffer, buffer, &buffer[total]);
		_currentOutBufSize = total;

		return total==0?EOF:this->buffer[0];
	}

protected:
	std::istream* in;
	unsigned char buffin[GZBUFSIZ];
	unsigned char buffout[GZBUFSIZ];
	char* buffer;
	z_stream strm;
	int status_flag;
	int _currentOutBufSize;
};



#endif /* INFLATESTREAMBUF_H_ */
