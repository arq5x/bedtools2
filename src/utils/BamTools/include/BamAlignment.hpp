#ifndef __HTSLIBPP_BAMALIGNEMNT_HPP__
#define __HTSLIBPP_BAMALIGNEMNT_HPP__
#include <stdint.h>
#include <memory>
#include <cstring>
namespace BamTools {
	static std::string _mkstr(const uint8_t* what) { return std::string((const char*)what + 1); }
	struct CigarOp {
	  
		char     Type;   //!< CIGAR operation type (MIDNSHPX=)
		uint32_t Length; //!< CIGAR operation length (number of bases)
		
		//! constructor
		CigarOp(const char type = '\0', 
				const uint32_t& length = 0)
			: Type(type)
			, Length(length) 
		{ }
	};
	const char cigar_ops_as_chars[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B' };

	class BamAlignment {

		struct _Bam {
			bam1_t* bam;
			~_Bam() { bam_destroy1(bam); }
			_Bam(bam1_t* b) : bam(b) {}
		};

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
				const uint8_t* data = bam_aux_get(this->_parent._bam->bam, tag.c_str());
				if(nullptr == data) return false;
				this->_dest = make_value(data);
				return true;
			}
		};

		std::shared_ptr<_Bam> _bam;
	public:

		const bam1_t* HtsObj() const 
		{
			return _bam->bam;
		}
		std::string Name;
		std::string Filename;
		std::vector<CigarOp> CigarData;

		void InitCigarData()
		{
			CigarData.clear();

			uint32_t *cigar = bam_get_cigar(_bam->bam);  
			uint32_t num_cigar_elements = _bam->bam->core.n_cigar;

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
			Qualities = "";
			static const char* base2char = "=ACMGRSVTWYHKDBN";
			const uint8_t* seq = bam_get_seq(_bam->bam);
			
			for(unsigned i = 0; i < QuerySequenceLength; i ++)
			{
				char cur_base = base2char[bam_seqi(seq, i)];
				QueryBases.push_back(cur_base);
			}

			const uint8_t* qual = bam_get_qual(_bam->bam);

			if(qual[0] == 0xffu)
				Qualities.resize(QuerySequenceLength, -1);
			else for(unsigned i = 0; i < QuerySequenceLength; i ++)
				Qualities.push_back((char)(33 + qual[i]));
			SupportData.AllCharData = std::string(_bam ? (char*)_bam->bam->data : "", _bam ? _bam->bam->l_data : 0);
		}

#include <BamAlignment.mapping.hpp>

		struct _SupportData {
			_QueryNameLength_t& QueryNameLength;
			_QuerySequenceLength_t& QuerySequenceLength;
			_NumCigarOperations_t& NumCigarOperations;
			uint32_t& BlockLength;
			std::string AllCharData; /* TODO(haohou): dealing with this dummy string, probably based on _ptr->data */
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
		

		BamAlignment() : _bam(NULL),  BlockLength(0), SupportData(*this){}

		BamAlignment(const std::string& filename, bam1_t* bam, uint32_t size = 0) : 
			_bam(new _Bam(bam)), 
			Name(bam_get_qname(bam)),
			Filename(filename),
			BlockLength(size),
			SupportData(*this)
		{
			setup(bam);
			InitCigarData();
			InitAdditionalData();
		}

		BamAlignment(const BamAlignment& ba) :
			_bam(ba._bam), 
			Filename(ba.Filename),
			SupportData(*this)
		{
			setup(ba);
			InitCigarData();
			InitAdditionalData();
		}

		const BamAlignment& operator = (const BamAlignment& ba)
		{
			_bam = ba._bam;
			Name = ba.Name;
			Filename = ba.Filename;
			setup(ba);
			BlockLength = ba.BlockLength;
			InitCigarData();
			InitAdditionalData();
			return *this;
		}

		void operator ()(const std::string filename, bam1_t* bam, uint32_t size = 0)
		{
			_bam = std::shared_ptr<_Bam>(new _Bam(bam));
			Filename = filename;
			Name = std::string(bam_get_qname(bam));
			setup(bam);
			BlockLength = size;
			InitCigarData();
			InitAdditionalData();
		}

		inline bool GetAlignmentFlag(uint32_t mask) const
		{
			return _bam->bam->core.flag & mask;
		}

		inline void SetAlignmentFlag(uint32_t mask, bool val)
		{
			if(val && !(_bam->bam->core.flag & mask))
				_bam->bam->core.flag |= mask;
			else if(val && (_bam->bam->core.flag & mask))
				_bam->bam->core.flag &= ~mask;
		}

#define FLAG_ACCESSOR(name, mask) \
		inline bool Is##name() const { return GetAlignmentFlag(mask); }; \
		inline void SetIs##name(bool val) { SetAlignmentFlag(mask, val); };

		inline bool IsMapped() const
		{
			return !(_bam->bam->core.flag & BAM_FUNMAP);
		}

		bool IsMateMapped() const
		{
			return !(_bam->bam->core.flag & BAM_FMUNMAP);
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

			return bam_endpos(_bam->bam);
		}

		template <typename T>
		bool GetTag(const std::string& tag, T& destination) const
		{

			_TagGetter<std::string, T> string_getter(*this, destination);
			if(string_getter)
				return string_getter(tag, _mkstr);
			return false;
		}

		bool HasTag(const std::string& tag) const
		{
			const uint8_t *aux_ptr = bam_aux_get(_bam->bam, tag.c_str());
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
				bam_aux_append(this->_bam->bam, tag.c_str(), 'Z', str->size() + 1, (uint8_t*)str->c_str());
				return true;
			}
			return false;
		}

	};
}
#endif
