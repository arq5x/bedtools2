#ifndef __HTSLIBPP_SAMHEADER_H__
#define __HTSLIBPP_SAMHEADER_H__
#include <stdint.h>
#include <string>
#include <htslib/sam.h>
#include <vector>
#include <string.h>
#include <api/BamAux.h>
#include <stdlib.h>


#ifdef WITH_HTS_CB_API
#include <htslib/hfile.h>
namespace htslib_future {
	static inline ssize_t buf_write(void* data, const void* buf, size_t sz) 
	{
		std::string* buffer = static_cast<std::string*>(data);

		buffer->append((const char*)buf, sz);

		return sz;
	}
	static inline int sam_hdr_rebuild(bam_hdr_t* hdr) 
	{
		std::string* buffer = new std::string();
		hFILE_callback_ops ops;
		memset(&ops, 0, sizeof(ops));
		ops.write = buf_write;
		ops.cb_data = buffer;
		samFile* fp = hts_open_callback(NULL, &ops, "w");

		sam_hdr_write(fp, hdr);

		hts_close(fp);

		hdr->l_text = (uint32_t)buffer->length();
		hdr->text = strdup(buffer->c_str());
		delete buffer;

		return 0;
	}
}
#endif

namespace BamTools {

	typedef std::vector<RefData> RefVector;
	
	static const std::string _defulat_sort_order = "coordinate";
	static const std::string _defualt_version = "1.4";
	static const std::string _default_group_order = "unknown";

	class SamHeader {
		bam_hdr_t* _header;
		std::string _filename;
	public:

		std::string SortOrder;
		std::string Version;
		std::string GroupOrder;

		SamHeader() : _header(NULL), _filename(""), SortOrder(_defulat_sort_order), Version(_defualt_version), GroupOrder(_default_group_order) {}

		SamHeader(const std::string& filename, bam_hdr_t* hdr) 
			: _header(hdr), _filename(filename), SortOrder(_defulat_sort_order), Version(_defualt_version),  GroupOrder(_default_group_order) 
		{
			char* buf = (char*)malloc(1024);
#ifdef WITH_HTS_CB_API
			htslib_future::sam_hdr_rebuild(_header);
#endif
			size_t len = _header->text == NULL ? 0 : strlen(_header->text);
			size_t sz = 0;
			for(size_t i = 0; i <= len; i ++)
			{
				char ch = _header->text[i];
				if(ch != ' ' &&  ch != '\n' && ch != '\t' && ch != '\n' && ch != 0)
				{
					if(sz >= 1024) buf = (char*)realloc(buf, len + 1);
					buf[sz++] = ch;
				}
				else 
				{
					buf[sz] = 0;
					if(sz > 3 &&  memcmp(buf, "VN:", 3) == 0)
					{
						std::string version_num(buf + 3);
						Version = version_num;
					}
					else if(sz > 3 && memcmp(buf, "SO:", 3) == 0)
					{
						std::string order(buf + 3);
						SortOrder = order;
					}
					sz = 0;
				}
			}
			free(buf);
		}

		std::string GetHeaderText() const 
		{
			return _header == NULL || _header->text == NULL ? "" : _header->text;
		}

		void ParseHeaderText(const std::string& text)
		{
			_header = sam_hdr_parse((int)text.length(), text.c_str());
			_header->text = strdup(text.c_str());
			_header->l_text = (uint32_t)text.length();
		}
		void destory() 
		{
			if(_header != NULL) bam_hdr_destroy(_header);
		}

		RefVector GetReferenceData() const
		{
			RefVector ret;
			for(int i = 0; i < _header->n_targets; i ++)
			{
				ret.push_back(RefData(_header->target_name[i], _header->target_len[i]));
			}
			return ret;
		}

		const std::string ToString() 
		{

			return _header != NULL && _header->text != NULL ? _header->text : "";
		}
		
		bool HasVersion() 
		{
			return true;
		}

		bam_hdr_t* GetHeaderStruct() 
		{
			return _header;
		}

		const char* Filename() 
		{
			return _filename.c_str();
		}
	};
}
#endif
