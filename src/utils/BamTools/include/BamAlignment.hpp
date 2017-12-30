#ifndef __HTSLIBPP_BAMALIGNEMNT_HPP__
#define __HTSLIBPP_BAMALIGNEMNT_HPP__
#include <stdint.h>
#include <Property.hpp>
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

		PropertyMapping<BamAlignment, int32_t> RefID;
		PropertyMapping<BamAlignment, int32_t> Position;
		PropertyMapping<BamAlignment, int32_t> MatePosition;
		PropertyMapping<BamAlignment, int32_t> MateRefID;
		//PropertyMapping1<BamAlignment, uint8_t, const uint16_t> MapQuality;
		
		uint16_t MapQuality;
		std::string Name;
		std::string Filename;

		struct _CigarData {
			operator std::vector<CigarOp>() const
			{
				uint32_t *cigar = bam_get_cigar(_parent._bam->bam);  
				uint32_t num_cigar_elements = _parent._bam->bam->core.n_cigar;
				std::vector<CigarOp> cd;

				for ( uint32_t i = 0; i < num_cigar_elements; ++i )
				{
					char type = cigar_ops_as_chars[static_cast<int>(bam_cigar_op(cigar[i]))];
					uint32_t len = bam_cigar_oplen(cigar[i]);
					cd.push_back(CigarOp(type, len));
				}
				return cd;
			}
			_CigarData(const BamAlignment& parent) : _parent(parent) {}
		private:
			const BamAlignment& _parent;
		} CigarData;
		
		std::string QueryBases, AlignedBases, Qualities;

		

		BamAlignment() : _bam(NULL),  CigarData(*this){}

		BamAlignment(const std::string& filename, bam1_t* bam) : 
			_bam(new _Bam(bam)), 
			RefID(_bam->bam->core.tid),
			Position(_bam->bam->core.pos),
			MatePosition(_bam->bam->core.mpos),
			MateRefID(_bam->bam->core.mtid),
			//MapQuality(_bam->bam->core.qual),
			Name(bam_get_qname(bam)),
			Filename(filename),
			CigarData(*this)
		{
			MapQuality = _bam->bam->core.qual;
		}

		BamAlignment(const BamAlignment& ba) :
			_bam(ba._bam), 
			RefID(_bam->bam->core.tid),
			Position(_bam->bam->core.pos),
			MatePosition(_bam->bam->core.mpos),
			MateRefID(_bam->bam->core.mtid),
			MapQuality(_bam->bam->core.qual),
			Name(ba.Name),
			Filename(ba.Filename),
			CigarData(*this)
		{
			MapQuality = _bam->bam->core.qual;
		}

		const BamAlignment& operator = (const BamAlignment& ba)
		{
			_bam = ba._bam;
			RefID(ba.RefID);
			Position(ba.Position);
			MateRefID(ba.MateRefID);
			MatePosition(ba.MatePosition);
			MapQuality = ba.MapQuality;
			Name = ba.Name;
			Filename = ba.Filename;
			return *this;
		}

		void operator ()(const std::string filename, bam1_t* bam)
		{
			_bam = std::shared_ptr<_Bam>(new _Bam(bam));
			RefID(_bam->bam->core.tid);
			Position(_bam->bam->core.pos);
			MatePosition(_bam->bam->core.mpos);
			MateRefID(_bam->bam->core.mtid);
			MapQuality = _bam->bam->core.qual;
			Filename = filename;
			Name = std::string(bam_get_qname(bam));
		}

		inline bool IsMapped() const
		{
			return !(_bam->bam->core.flag & BAM_FUNMAP);
		}

		inline bool IsFirstMate() const
		{
			return _bam->bam->core.flag & BAM_FREAD1;
		}

		inline bool IsSecondMate() const
		{
			return _bam->bam->core.flag & BAM_FREAD2;
		}

		inline bool IsReverseStrand() const
		{
			return _bam->bam->core.flag & BAM_FREVERSE;
		}

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

		bool IsMateMapped() const
		{
			return !(_bam->bam->core.flag & BAM_FMUNMAP);
		}

		template <typename T>
		bool AddTag(const std::string& tag, const std::string& type, T& data)
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
