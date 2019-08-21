#ifndef __HTSLIBPP_SAMHEADER_H__
#define __HTSLIBPP_SAMHEADER_H__
#include <stdint.h>
#include <string>
#include <htslib/sam.h>
#include <vector>
#include <string.h>
#include <api/BamAux.h>
#include <stdlib.h>
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
			_header = bam_hdr_init();
			int sq = 0;
			const char* ptr = text.c_str();
			std::vector<uint32_t> target_size;
			std::vector<std::string> target_name;
			size_t sz = 0;
			char* buf = (char*)malloc(1024);
			for(;;ptr++)
			{
				if(*ptr &&  *ptr != ' ' &&  *ptr != '\n' && *ptr != '\t' && *ptr != '\n')
				{
					if(sz >= 1024) buf = (char*)realloc(buf, text.length() + 1);
					buf[sz++] = *ptr;
				}
				else
				{
					if(sz > 0)
					{
						buf[sz] = 0;
						if(sq == 1)
						{
							target_name.push_back(buf + 3);
							sq = 2;
						}
						else if(sq == 2)
						{
							target_size.push_back(atoi(buf + 3));
							sq = 0;
						}
						if(strcmp(buf, "@SQ") == 0)
							sq = 1;
					}
					sz = 0;
				}
				if(*ptr == 0) break;
			}

			free(buf);

			_header->n_targets = target_size.size();
			_header->l_text = text.length();
			_header->target_len = (uint32_t*)malloc(sizeof(uint32_t) * _header->n_targets);
			_header->target_name = (char**)malloc(sizeof(char*) * _header->n_targets);

			for(int i = 0; i < _header->n_targets; i++)
			{
				_header->target_len[i] = target_size[i];
				_header->target_name[i] = strdup(target_name[i].c_str());
			}

			_header->text = strdup(text.c_str());

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
