#ifndef __HTSLIBPP_BAMREADER_HPP__
#define __HTSLIBPP_BAMREADER_HPP__
#include <BamAlignment.hpp>
#include <Property.hpp>
#include <SamHeader.hpp>
#include <htslib/sam.h>
#include <stdint.h>
#include <string>
#include <queue>

namespace BamTools {
	class BamReader {
		struct _SamFile {
			_SamFile(samFile* fp, uint32_t _idx) : fp(fp), idx(_idx){}
			void destory() 
			{
				if(nullptr != fp) sam_close(fp);
			}
			samFile* fp;
			uint32_t idx;
		};

		struct _Comp {
			bool operator()(const std::pair<_SamFile, bam1_t*>& left, const std::pair<_SamFile, bam1_t*>& right) 
			{
				if(left.second->core.tid == right.second->core.tid)
				{
					return left.second->core.pos >= right.second->core.pos;
				}
				return left.second->core.tid >= right.second->core.tid;
			}
		};

		std::vector<_SamFile> _files;
		std::vector<SamHeader> _hdrs;
		std::priority_queue<std::pair<_SamFile, bam1_t*>, std::vector<std::pair<_SamFile, bam1_t*> >, _Comp> _queue;

		bool _read_sam_file(_SamFile file)
		{
			bam1_t *rec_ptr = bam_init1();
			if (sam_read1(file.fp, _hdrs[file.idx].GetHeaderStruct(), rec_ptr) >= 0)
			{
				_queue.push(std::make_pair(file, rec_ptr));
				return true;
			}
			bam_destroy1(rec_ptr);
			return false;
		}

	public:
		bool Open(const std::vector<std::string>& filenames) 
		{
			uint32_t i = 0;
			for(const auto& filename : filenames)
			{
				samFile* fp = sam_open(filename.empty() ? "-" : filename.c_str(), "rb");
				if(nullptr == fp) return false;
				bam_hdr_t* hdr = sam_hdr_read(fp);
				if(nullptr == hdr)
					return false;

				_files.push_back(_SamFile(fp, i ++));
				_hdrs.push_back(SamHeader(filename, hdr));

				_read_sam_file(_files[_files.size() - 1]);
			}

			return true;
		}

		bool Open(const std::string& filename) 
		{
			std::vector<std::string> vec;
			vec.push_back(filename);
			return Open(vec);
		}

		bool IsOpen() const 
		{
			return true;
		}

		const RefVector GetReferenceData() const
		{
			return _hdrs.at(0).GetReferenceData();
		}

		SamHeader GetHeader(int idx = 0) const
		{
			return _hdrs.at(idx);
		}

		std::string GetHeaderText(int idx = -1) const 
		{
			return _hdrs.at(idx).GetHeaderText();
		}
		std::string GetErrorString() const
		{
			return "FIXME: error string is not supported";
		}

		bool GetNextAlignment(BamAlignment& alignment)
		{
			if(_queue.empty()) return false;

			auto& top = _queue.top();

			_SamFile fp = top.first;
			
			alignment(_hdrs[fp.idx].Filename(), top.second);

			_queue.pop();

			_read_sam_file(fp);

			return true;
		}

		bool GetNextAlignmentCore(BamAlignment& alignment)
		{
			return GetNextAlignment(alignment);
		}

		void Close(void)
		{
			for(auto& sam : _files)
			{
				sam.destory();
			}

			for(auto& hdr : _hdrs)
			{
				hdr.destory();
			}

			while(!_queue.empty())
			{
				auto& item = _queue.top();
				bam_destroy1(item.second);
				_queue.pop();
			}
		}
	};
}
#endif
