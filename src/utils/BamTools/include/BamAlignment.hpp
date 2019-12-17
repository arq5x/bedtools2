#ifndef __HTSLIBPP_BAMALIGNEMNT_HPP__
#define __HTSLIBPP_BAMALIGNEMNT_HPP__
#include <stdint.h>
#include <memory>
#include <cstring>
#include <stdlib.h>
#include <htslib/hts_endian.h>
namespace BamTools {
	const static char cigar_ops_as_chars[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B' };
	static inline std::string _mkstr(const uint8_t* what) { return std::string((const char*)what + 1); }
	template <class T>
	static inline T _mk_numeric(const uint8_t* what) {
		switch(what[0]) {
			case 'd': 
				return le_to_double(what + 1);
			case 'f':
				return le_to_float(what + 1);
			case 'A':
			case 'C':
				return what[1];
			case 'c':
				return le_to_i8(what + 1);
			case 's':
				return le_to_i16(what + 1);
			case 'S':
				return le_to_u16(what + 1);
			case 'i':
				return le_to_i32(what + 1);
			case 'I':
				return le_to_u32(what + 1);
			default:
				return 0;
		}
	}
	struct CigarOp {
	  
		char     Type;   //!< CIGAR operation type (MIDNSHPX=)
		uint32_t Length; //!< CIGAR operation length (number of bases)
		
		//! constructor
		CigarOp(const char type = '\0', 
				const uint32_t& length = 0)
			: Type(type)
			, Length(length) 
		{ }

		uint32_t to_htslib_repr() const
		{
			static bool init_flag = false;
			static int cigar_code[256] = {};

			if(!init_flag) 
			{
				init_flag = true;
				memset(cigar_code, -1, sizeof(cigar_code));
				for(unsigned i = 0; i < sizeof(cigar_ops_as_chars)/sizeof(*cigar_ops_as_chars); i++)
					cigar_code[(size_t)cigar_ops_as_chars[i]] = (int)i;
			}

			int code = cigar_code[(size_t)Type];
			if(code < 0) return 0;

			return (Length << 4) + code;
		}
	};

	class BamAlignment {

		bam1_t _bam;

		template <typename Destination>
		struct _TagGetterBase 
		{
			_TagGetterBase(const BamAlignment& parent, Destination& dest) : _dest(dest), _parent(parent) {}
			Destination& _dest;
			const BamAlignment& _parent;
		};

		template <typename Left, typename Right> 
		struct _TagGetter : _TagGetterBase<Right>
		{
			operator bool() { return false; }

			_TagGetter(const BamAlignment& parent, Right& right) : _TagGetterBase<Right>(parent, right) {}

			template <class MakeValue>
			bool operator()(const std::string tag, MakeValue make_value)
			{
				return false;
			}
		};
		
		template <typename Left>
		struct _TagGetter<Left, Left> : _TagGetterBase<Left>
		{
			operator bool() { return true; }
			
			_TagGetter(const BamAlignment& parent, Left& right) : _TagGetterBase<Left>(parent, right){};
			
			template <class MakeValue>
			bool operator ()(const std::string& tag, MakeValue make_value)
			{
				const uint8_t* data = bam_aux_get(this->_parent.HtsObj(), tag.c_str());
				if(nullptr == data) return false;
				this->_dest = make_value(data);
				return true;
			}
		};

		static inline void* _ensure_data_chunk(bam1_t* bam, size_t ofs, size_t old_size, size_t new_size) 
		{
			if(ofs + old_size > (size_t)bam->l_data) return NULL;

			if(bam->data == NULL || bam->m_data < bam->l_data + new_size - old_size)
			{
				size_t next_size = (bam->l_data + new_size - old_size + 31) & ~(size_t)31;
				if(NULL == (bam->data = (uint8_t*)realloc(bam->data, next_size))) return NULL;
				bam->m_data = next_size;
			}

			memmove(bam->data + ofs + old_size, bam->data + ofs + new_size, bam->l_data - ofs - old_size);

			bam->l_data = bam->l_data + new_size - old_size;

			return bam->data + ofs;
		}

	public:

		const bam1_t* HtsObj() const 
		{
			return &_bam;
		}

		bam1_t* HtsObj2() 
		{
			return &_bam;
		}

		std::string Name;
		std::string Filename;
		std::vector<CigarOp> CigarData;

		void InitCigarData()
		{
			CigarData.clear();

			uint32_t *cigar = bam_get_cigar(&_bam);  
			uint32_t num_cigar_elements = _bam.core.n_cigar;

			for ( uint32_t i = 0; i < num_cigar_elements; ++i )
			{
				char type = cigar_ops_as_chars[static_cast<int>(bam_cigar_op(cigar[i]))];
				uint32_t len = bam_cigar_oplen(cigar[i]);
				CigarData.push_back(CigarOp(type, len));
			}
		}
		
		/* TODO: fix all this dummy strings What is TagData? Not from file? */
		std::string AlignedBases, Qualities, ErrorString, TagData;
		uint32_t BlockLength;

		std::string QueryBases;

		void InitAdditionalData()
		{
			QueryBases = "";
			static const char* base2char = "=ACMGRSVTWYHKDBN";
			const uint8_t* seq = bam_get_seq(&_bam);
			
			for(unsigned i = 0; i < QuerySequenceLength; i ++)
			{
				char cur_base = base2char[bam_seqi(seq, i)];
				QueryBases.push_back(cur_base);
			}

			Qualities = "";
			
			const uint8_t* qual = bam_get_qual(&_bam);

			if(_bam.core.l_qseq == 0 || qual[0] == 0xffu)
				Qualities.resize(QuerySequenceLength, 33);
			else for(unsigned i = 0; i < QuerySequenceLength; i ++)
				Qualities.push_back((char)(33 + qual[i]));
		}

		bool SyncExtraData() const
		{
#define BAM_DATA_OFFSET(what) ((size_t)(((uint8_t*)bam_get_##what(&_bam)) - ((uint8_t*)_bam.data)))
			void* qname_buf = _ensure_data_chunk((bam1_t*)&_bam, BAM_DATA_OFFSET(qname), _bam.core.l_qname, Name.size() + 1);
			if(NULL == qname_buf) return false;
			memcpy(qname_buf, Name.c_str(), Name.size() + 1);
			((bam1_t*)&_bam)->core.l_qname = Name.size() + 1;

			uint32_t* cigar_buf = (uint32_t*)_ensure_data_chunk((bam1_t*)&_bam, BAM_DATA_OFFSET(cigar), _bam.core.n_cigar * sizeof(uint32_t), sizeof(uint32_t) * CigarData.size());
			if(NULL == cigar_buf) return false;
			for(auto& cigar_op : CigarData) 
				*(cigar_buf++) = cigar_op.to_htslib_repr();
			((bam1_t*)&_bam)->core.n_cigar = CigarData.size();

			// TODO: also sync the quality and qseq
#undef BAM_DATA_OFFSET

			return true;
		}

#include <BamAlignment.mapping.hpp>

		struct _SupportData {
			_QueryNameLength_t& QueryNameLength;
			_QuerySequenceLength_t& QuerySequenceLength;
			_NumCigarOperations_t& NumCigarOperations;
			uint32_t& BlockLength;
			struct {
				const char* begin;
				size_t      size;
				const char* c_str() const { return begin; }
				size_t      length() const { return size; }
				void clear() { begin = NULL; size = 0; }
			} AllCharData;
			bool HasCoreOnly;  /* TODO(haohou): populate the string data */
			/* TODO(haohou): Tag2Cigar?  Real Cigar ? */
			_SupportData(BamAlignment& parent): 
				QueryNameLength(parent.QueryNameLength),
				QuerySequenceLength(parent.QuerySequenceLength),
				NumCigarOperations(parent.NumCigarOperations),
				BlockLength(parent.BlockLength),
				HasCoreOnly(false)
			{
			}
		} SupportData;
	
		~BamAlignment() 
		{
			free(_bam.data);
		}

		BamAlignment() : BlockLength(0), SupportData(*this)
		{
			memset(&_bam, 0, sizeof(_bam));
		}


		BamAlignment(const std::string& filename, const bam1_t* bam, uint32_t size = 0) : 
			Name(bam_get_qname(bam)),
			Filename(filename),
			BlockLength(size),
			SupportData(*this)
		{
			memset(&_bam, 0, sizeof(_bam));
			bam_copy1(&_bam, bam);
			InitCigarData();
			InitAdditionalData();
		}

		BamAlignment(const BamAlignment& ba) :
			Filename(ba.Filename),
			SupportData(*this)
		{
			memset(&_bam, 0, sizeof(_bam));
			bam_copy1(&_bam, ba.HtsObj());
			Name = ba.Name;
			CigarData = ba.CigarData;
			QueryBases = ba.QueryBases;
			Qualities = ba.Qualities;
			SupportData.AllCharData = ba.SupportData.AllCharData;
		}

		const BamAlignment& operator = (const BamAlignment& ba)
		{
			Name = ba.Name;
			Filename = ba.Filename;
			bam_copy1(&_bam, ba.HtsObj());
			BlockLength = ba.BlockLength;
			CigarData = ba.CigarData;
			QueryBases = ba.QueryBases;
			Qualities = ba.Qualities;
			SupportData.AllCharData = ba.SupportData.AllCharData;
			return *this;
		}

		void operator ()(const std::string filename, const bam1_t* bam, uint32_t size = 0, bool copy = true)
		{
			Filename = filename;
			Name = std::string(bam_get_qname(bam));
			BlockLength = size;
			if(copy)
				bam_copy1(&_bam, bam);
			else
			{
				free(_bam.data);
				_bam = *bam;
			}

			SupportData.AllCharData.begin = (const char*)_bam.data;
			SupportData.AllCharData.size = _bam.l_data;

			InitCigarData();
		}

		inline bool GetAlignmentFlag(uint32_t mask) const
		{
			return _bam.core.flag & mask;
		}

		inline void SetAlignmentFlag(uint32_t mask, bool val)
		{
			if(val && !(_bam.core.flag & mask))
				_bam.core.flag |= mask;
			else if(val && (_bam.core.flag & mask))
				_bam.core.flag &= ~mask;
		}

#define FLAG_ACCESSOR(name, mask) \
		inline bool Is##name() const { return GetAlignmentFlag(mask); }; \
		inline void SetIs##name(bool val) { SetAlignmentFlag(mask, val); };

		inline bool IsMapped() const
		{
			return !(_bam.core.flag & BAM_FUNMAP);
		}

		inline bool IsMateMapped() const
		{
			return !(_bam.core.flag & BAM_FMUNMAP);
		}


		FLAG_ACCESSOR(FirstMate, BAM_FREAD1);
		FLAG_ACCESSOR(SecondMate, BAM_FREAD2);
		FLAG_ACCESSOR(ReverseStrand, BAM_FREVERSE);
		FLAG_ACCESSOR(MateReverseStrand, BAM_FMREVERSE);
		FLAG_ACCESSOR(ProperPair, BAM_FPROPER_PAIR);
		FLAG_ACCESSOR(Paired, BAM_FPAIRED);
		FLAG_ACCESSOR(Duplicate, BAM_FDUP);
		FLAG_ACCESSOR(FailedQC, BAM_FQCFAIL);

		int GetEndPosition(bool usePadded = false, bool closedInterval = false) const
		{
			if(usePadded || closedInterval) 
			{
				fprintf(stderr, "FIXME: we current do not support the usePadded and closedInterval param");
			}

			return bam_endpos(&_bam);
		}

		template <typename T>
		bool GetTag(const std::string& tag, T& destination) const
		{

			{
				_TagGetter<std::string, T> string_getter(*this, destination);
				if(string_getter)
					return string_getter(tag, _mkstr);
			}

#define _GET_NUM(type) \
			{\
				_TagGetter<type, T> getter(*this, destination);\
				if(getter)\
					return getter(tag, _mk_numeric<type>);\
			}
			_GET_NUM(uint8_t);
			_GET_NUM(int8_t);
			_GET_NUM(uint16_t);
			_GET_NUM(int16_t);
			_GET_NUM(uint32_t);
			_GET_NUM(int32_t);
			_GET_NUM(double);
			_GET_NUM(float);

			return false;
		}

		bool HasTag(const std::string& tag) const
		{
			const uint8_t *aux_ptr = bam_aux_get(&_bam, tag.c_str());
			if ( aux_ptr == nullptr )
				return false;
			return true;
		}
		
		template <typename T>
		bool AddTag(const std::string& tag, const std::string& type, T data)
		{
			_TagGetter<std::string, T> string_getter(*this, data);
			if(string_getter && "Z" == type)
			{
				std::string* str = static_cast<std::string*>(&data);
				bam_aux_append(&_bam, tag.c_str(), 'Z', str->size() + 1, (uint8_t*)str->c_str());
				return true;
			}
			return false;
		}

	};
}
#endif
