#ifndef __HTSLIBPP_BAM_WRITER_HPP__
#define __HTSLIBPP_BAM_WRITER_HPP__
#include <string>
#include <htslib/sam.h>
#include <SamHeader.hpp>
#include <BamAlignment.hpp>
namespace BamTools {
	class BamWriter {
		samFile* _fp;
		SamHeader _hdr;
	public:
		bool Open(const std::string& filename, const std::string& samHeaderText, const RefVector& referenceSequences, bool binary = true)
		{
			_fp = sam_open(filename.empty() ? "-" : filename.c_str(), binary ? "wb" : "w");
			if(_fp == nullptr) return false;

			_hdr.ParseHeaderText(samHeaderText);
			sam_hdr_write(_fp, _hdr.GetHeaderStruct());
			
			return true;
		}
		void Close(void)
		{
			if(nullptr != _fp) sam_close(_fp);
		}
		void SaveAlignment(const BamAlignment& al)
		{
			const bam1_t* bam = al.HtsObj();
			sam_write1(_fp, _hdr.GetHeaderStruct(), bam);
		}
	};
}
#endif
