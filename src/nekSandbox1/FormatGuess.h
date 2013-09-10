///*
// * FormatGuess.h
// *
// *  Created on: Mar 11, 2013
// *      Author: nek3d
// */
//
//#ifndef FORMATGUESS_H_
//#define FORMATGUESS_H_
//
//#include "PushBackStream.h"
//
//class FormatGuess
//{
//	std::string name;
//public:
//	FormatGuess(const char* name)
//	:name(name)
//	{}
//	const char* format() const
//	{
//		return name.c_str();
//	}
//	virtual bool guess(PushBackStreamBuf* buf)=0;
//};
//
//class VCFGuess:public FormatGuess
//{
//public:
//	VCFGuess()
//	:FormatGuess("VCF")
//	{}
//
//	virtual bool guess(PushBackStreamBuf* buf)
//	{
//		const std::string fileformat("##fileformat=");
//		bool is_vcf=true;
//		std::ostringstream os;
//		for(std::size_t i=0;i< fileformat.size();++i)
//		{
//			int c=buf->sbumpc();
//			if(c!=-1) os << (char)c;
//			if(c!=fileformat[i]) {is_vcf=false;break;}
//		}
//		buf->push_back(os.str());
//		return is_vcf;
//	}
//};
//
//
//class XMLGuess:public FormatGuess
//{
//public:
//
//	XMLGuess()
//	:FormatGuess("XML")
//	{}
//
//	virtual bool guess(PushBackStreamBuf* buf)
//	{
//		const std::string xmlheader("<?xml ");
//		bool is_xml=true;
//		std::ostringstream os;
//		for(std::size_t i=0;i< xmlheader.size();++i)
//		{
//			int c=buf->sbumpc();
//			if(c!=-1) os << (char)c;
//			if(c!=xmlheader[i]) {is_xml=false;break;}
//		}
//		buf->push_back(os.str());
//		return is_xml;
//	}
//};
//
//
//#endif /* FORMATGUESS_H_ */
