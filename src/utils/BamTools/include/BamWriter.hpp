#ifndef __HTSLIBPP_BAM_WRITER_HPP__
#define __HTSLIBPP_BAM_WRITER_HPP__
#include <string>
#include <stdexcept>
#include <htslib/sam.h>
#include <SamHeader.hpp>
#include <BamAlignment.hpp>

extern const char* cram_reference;

namespace BamTools {
	class BamWriter {
		samFile* _fp;
		SamHeader _hdr;
	public:
		enum CompressionMode {
			Compressed,
			Uncompressed
		};

		bool Open(const std::string& filename, const std::string& samHeaderText, const RefVector& referenceSequences, refs_t* reference)
		{
			const char* ref_file = cram_reference ? cram_reference : getenv("CRAM_REFERENCE");

			_fp = sam_open(filename.empty() || filename == "stdout" ? "-" : filename.c_str(), ref_file ? "wc" : "wb");
			if(_fp == nullptr) return false;

			if(hts_set_opt(_fp, CRAM_OPT_SHARED_REF, reference) < 0)
				return false;

			_hdr.ParseHeaderText(samHeaderText);

			return sam_hdr_write(_fp, _hdr.GetHeaderStruct()) >= 0;
		}

		bool Open(const std::string& filename, const std::string& samHeaderText, const RefVector& referenceSequences, bool binary = true)
		{
			_fp = sam_open(filename.empty() || filename == "stdout"? "-" : filename.c_str(), binary ? "wb" : "w");
			if(_fp == nullptr) return false;

			_hdr.ParseHeaderText(samHeaderText);
			return sam_hdr_write(_fp, _hdr.GetHeaderStruct()) >= 0;
		}

		void Close(void)
		{
			if(nullptr != _fp) sam_close(_fp);
		}

		void SaveAlignment(const BamAlignment& al, bool sync_data = false)
		{
			if(sync_data)
				al.SyncExtraData();
			const bam1_t* bam = al.HtsObj();
			if (sam_write1(_fp, _hdr.GetHeaderStruct(), bam) < 0)
				throw std::runtime_error("can't write alignment record");
		}
		void SetCompressionMode(CompressionMode mode)
		{
			/* TODO(haohou) implemnt this */
		}
	};
}
#endif
