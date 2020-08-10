#ifndef __HTSLIBPP_BAMREADER_HPP__
#define __HTSLIBPP_BAMREADER_HPP__
#include <BamAlignment.hpp>
#include <SamHeader.hpp>
#include <htslib/cram.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <stdint.h>
#include <string>
#include <queue>
#include <istream>

extern const char* cram_reference;
namespace BamTools {
#ifdef WITH_HTS_CB_API
	struct stream_data_t {
		std::istream& stream;
		stream_data_t(std::istream& is) : stream(is){}
		
		static ssize_t read(void* data, void* resbuf, size_t sz)
		{
			stream_data_t *stream_data = static_cast<stream_data_t*>(data);
			stream_data->stream.read((char*)resbuf, sz);
			if(!stream_data->stream)
				sz = (ssize_t)stream_data->stream.gcount();
			return (ssize_t)sz;
		}
		
		static int close(void* data)
		{
			delete static_cast<stream_data_t*>(data);
			return 0;
		}
	};



#endif
	class BamReader {
		struct _SamFile {
			_SamFile(samFile* fp, uint32_t _idx, BamReader* reader) : 
				fp(fp), idx(_idx), ip(NULL), it(NULL), 
				reader(reader), has_range(false)
			{}
			~_SamFile() 
			{
				if(nullptr != ip) hts_idx_destroy(ip);
				if(nullptr != fp) sam_close(fp);
				if(nullptr != it) hts_itr_destroy(it);
			}
			
			bool load_index(const char* filename) 
			{
				if(ip == NULL && NULL == (ip = sam_index_load(fp, filename)))
					return false;
				return true;
			}

			bool refresh_idx()
			{
				if(NULL != it) 
					hts_itr_destroy(it);
				it = NULL;
				if(tid_l > tid_r || (tid_l == tid_r && ofs_l > ofs_r)) return false;
				if(NULL == (it = sam_itr_queryi(ip, tid_l, ofs_l, tid_r != tid_l ? reader->GetReferenceData().at(tid_l).RefLength : ofs_r + 1)))
					return false;
				return true;
			}
			bool set_range(BamRegion& reg)
			{
				tid_l = reg.LeftRefID;
				tid_r = reg.RightRefID;
				ofs_l = reg.LeftPosition;
				ofs_r = reg.RightPosition;
				has_range = true;
				return refresh_idx();
			}

			bool update_range()
			{
				tid_l ++;
				ofs_l = 0;
				return refresh_idx();
			}

			refs_t* get_reference()
			{
				return cram_get_refs(fp);
			}

			void set_cram_reqd_fields(int reqd_fields_code)
			{
				// tell HTSLIB what fields to read via
				// specifying reqd_fields_code
				hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, reqd_fields_code);
			}

			samFile* fp;
			uint32_t idx;
			hts_idx_t* ip;
			uint32_t cram_fields_code;
			
			int tid_l, tid_r;
			int ofs_l, ofs_r;
			hts_itr_t* it;

			BamReader* reader;
			bool has_range;
			bam1_t buffer;
		};

		struct _MetaData {
			_SamFile* file;
			uint32_t size;
			_MetaData(_SamFile* _file, uint32_t size) : file(_file)
			{
				this->size = size;
			}
		};


		struct _Comp {
			bool operator()(const std::pair<_MetaData, bam1_t*>& left, const std::pair<_MetaData, bam1_t*>& right) 
			{
				if(left.second->core.tid == right.second->core.tid)
				{
					return left.second->core.pos >= right.second->core.pos;
				}
				return left.second->core.tid >= right.second->core.tid;
			}
		};

		std::vector<_SamFile*> _files;
		std::vector<SamHeader> _hdrs;
		std::priority_queue<std::pair<_MetaData, bam1_t*>, std::vector<std::pair<_MetaData, bam1_t*> >, _Comp> _queue;
		std::string _error_str;

		bool _read_sam_file(_SamFile* file)
		{
			memset(&file->buffer, 0, sizeof(file->buffer));
			bam1_t* rec_ptr = &file->buffer;

			int read_rc;
			if(file->has_range)
			{
				read_rc = file->it != NULL ? sam_itr_next(file->fp, file->it, rec_ptr) : -1;
				if(file->it != NULL && read_rc == -1) 
				{
					file->update_range();
					read_rc = file->it != NULL ? sam_itr_next(file->fp, file->it, rec_ptr) : -1;
				}
			}
			else 
				read_rc = sam_read1(file->fp, _hdrs[file->idx].GetHeaderStruct(), rec_ptr);

			_error_str = "";

			if (read_rc >= 0)
			{
				_queue.push(std::make_pair(_MetaData(file, (uint32_t)(read_rc - 4)) , rec_ptr));
				return true;
			}
			else if(read_rc < -1)
			{
				_error_str = "Htslib error"; 
			}

			return false;
		}

		bool _Open_impl(uint32_t idx, samFile* fp, const std::string& filename = "-")
		{
			if(nullptr == fp) return false;
			if(fp->format.format != htsExactFormat::bam && fp->format.format != htsExactFormat::cram && fp->format.format != htsExactFormat::sam)
				return false;

			_SamFile* sam_file = new _SamFile(fp, idx, this);

			if(fp->format.format == htsExactFormat::cram)
			{
				const char* ref_file = cram_reference ? cram_reference : getenv("CRAM_REFERENCE");
				if(NULL == ref_file || hts_set_fai_filename(fp, ref_file) == -1)
				{
					// If we are not able to load a reference 
					// we will punt and just load the core alignment data.
					sam_file->set_cram_reqd_fields(2559);
				}
			}
			
			bam_hdr_t* hdr = sam_hdr_read(fp);
			if(nullptr == hdr)
				return false;


			_files.push_back(sam_file);
			_hdrs.push_back(SamHeader(filename, hdr));

			if(_read_sam_file(_files[_files.size() - 1]) || GetErrorString() == "")
				return true;

			Close();

			return false;
		}

		bool _Open_impl(uint32_t idx, const std::string& filename)
		{
			samFile* fp = NULL;
			fp = sam_open(filename.empty() || filename == "stdin" ? "-" : filename.c_str(), "rb");
			return _Open_impl(idx, fp, filename);
		}

	public:
		bool Open(const std::vector<std::string>& filenames) 
		{
			uint32_t i = 0;
			for(const auto& filename : filenames)
				if(_Open_impl(i, filename))
					i ++;
				else
					return false;
			return true;
		}

		bool Open(const std::string& filename) 
		{
			return _Open_impl(0, filename);
		}

#ifdef WITH_HTS_CB_API
		bool OpenStream(std::istream* is)
		{
			if(nullptr == is) return false;
			stream_data_t* cb_data = new stream_data_t(*is);
			hFILE_callback_ops ops;
			memset(&ops, 0, sizeof(ops));
			ops.read = stream_data_t::read;
			ops.close = stream_data_t::close;
			ops.cb_data = cb_data;
			samFile* fp = hts_open_callback(NULL, &ops, "rb");
			return _Open_impl(0, fp);
		}
#endif

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

		std::string GetHeaderText(int idx = 0) const 
		{
			return _hdrs.at(idx).GetHeaderText();
		}
		std::string GetErrorString() const
		{
			return _error_str;
		}

		bool GetNextAlignment(BamAlignment& alignment)
		{
			if(GetNextAlignmentCore(alignment))
			{
				alignment.InitAdditionalData();
				return true;
			}
			return false;
		}

		bool GetNextAlignmentCore(BamAlignment& alignment)
		{
			if(_queue.empty()) return false;

			auto& top = _queue.top();

			_SamFile* fp = top.first.file;
			
			alignment(_hdrs[fp->idx].Filename(), top.second, top.first.size, false);

			_queue.pop();

			_read_sam_file(fp);

			return true;
		}

		void Close(void)
		{
			for(auto& sam : _files)
			{
				delete sam;
			}

			for(auto& hdr : _hdrs)
			{
				hdr.destory();
			}

			while(!_queue.empty())
			{
				auto& item = _queue.top();
				free(item.second->data);
				_queue.pop();
			}
		}

		int GetReferenceID(const std::string& refname)
		{
			if(_hdrs.size() == 0) return -1;

			bam_hdr_t* bh = _hdrs[0].GetHeaderStruct();

			for(int i = 0; i < bh->n_targets; i ++)
				if(strcmp(refname.c_str(), bh->target_name[i]) == 0)
					return i;
			return -1;
		}

		void LocateIndexes() 
		{
			for(auto& sam: _files)
			{
				if(!sam->load_index(_hdrs[sam->idx].Filename()))
				{
					/* TODO(haohou): Load failure */
				}
			}
		}

		bool HasIndexes()
		{
			for(auto& sam: _files)
			{
				if(NULL == sam->ip)
					return false;
			}
			return true;
		}

		bool SetRegion(BamRegion& region)
		{
			for(auto& sam: _files)
			{
				if(!sam->set_range(region))
					return false;
			}

			_queue = std::priority_queue<std::pair<_MetaData, bam1_t*>, std::vector<std::pair<_MetaData, bam1_t*> >, _Comp>();

			for(auto& sam: _files)
			{
				_error_str = "";
				if(!_read_sam_file(sam) && _error_str != "")
					return false;
			}

			return true;
		}

		void SetCramReqdFieldsCode(int ReqdFieldsCode)
		{
			// tell HTSLIB what fields to read via
			// specifying reqd_fields_code
			for(auto& file: _files)
			{
				file->set_cram_reqd_fields(ReqdFieldsCode);
			}
		}

		refs_t* GetReference(int idx = 0)
		{
			return _files.at(idx)->get_reference();
		}

	};
}
#endif
