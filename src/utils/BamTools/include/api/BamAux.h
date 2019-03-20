#include <string>

#ifndef BAMAUX_H
#define BAMAUX_H
namespace BamTools {
	struct RefData
	{

		std::string RefName;  //!< name of reference sequence
		int32_t RefLength;    //!< length of reference sequence

		//! constructor
		RefData(const std::string& name = std::string(), const int32_t& length = 0)
			: RefName(name)
			, RefLength(length)
		{}
	};
	struct BamRegion
	{

		int LeftRefID;      //!< reference ID for region's left boundary
		int LeftPosition;   //!< position for region's left boundary
		int RightRefID;     //!< reference ID for region's right boundary
		int RightPosition;  //!< position for region's right boundary

		//! constructor
		BamRegion(const int& leftID = -1, const int& leftPos = -1, const int& rightID = -1,
				  const int& rightPos = -1)
			: LeftRefID(leftID)
			, LeftPosition(leftPos)
			, RightRefID(rightID)
			, RightPosition(rightPos)
		{}

		//! copy constructor
		BamRegion(const BamRegion& other)
			: LeftRefID(other.LeftRefID)
			, LeftPosition(other.LeftPosition)
			, RightRefID(other.RightRefID)
			, RightPosition(other.RightPosition)
		{}

		//! Clears region boundaries
		void clear()
		{
			LeftRefID = -1;
			LeftPosition = -1;
			RightRefID = -1;
			RightPosition = -1;
		}

		//! Returns true if region has a left boundary
		bool isLeftBoundSpecified() const
		{
			return (LeftRefID >= 0 && LeftPosition >= 0);
		}

		//! Returns true if region boundaries are not defined
		bool isNull() const
		{
			return (!isLeftBoundSpecified() && !isRightBoundSpecified());
		}

		//! Returns true if region has a right boundary
		bool isRightBoundSpecified() const
		{
			return (RightRefID >= 0 && RightPosition >= 1);
		}
	};
}
#endif
