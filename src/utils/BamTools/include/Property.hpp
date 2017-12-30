#ifndef __HTSLIBPP_PROPERTY_HPP__
#define __HTSLIBPP_PROPERTY_HPP__
#include <string>
namespace BamTools {
	template <class Host, class Type> 
	struct PropertyMapping {
		operator Type&() const
		{
			return *_target;
		}
		const Type& operator = (const Type& order) const
		{
			return NULL != _target ? (*_target = order) : order;
		}
		const Type& operator = (const PropertyMapping that) const
		{
			return *this = *that._target;
		}
		void operator ()(PropertyMapping& that)
		{
			_target = that._target;
		}
		void operator ()(Type& obj) 
		{
			_target = &obj;
		}
		PropertyMapping(Type& obj) : _target(&obj) {}
		PropertyMapping() : _target(NULL) {}
	private:
		Type* _target;
	};

	template <class Host, class SrcT, class DestT>
	struct PropertyMapping1 {
		operator DestT&() const 
		{
			return reinterpret_cast<DestT>(*_target);
		}
		const DestT& operator = (const SrcT& order) const 
		{
			return NULL != _target ? (*_target = reinterpret_cast<SrcT>(order)) : order;
		}
		void operator()(PropertyMapping1& that)
		{
			_target = that._target;
		}
		void operator ()(SrcT& obj)
		{
			_target = &obj;
		}
		PropertyMapping1(SrcT& obj) : _target(&obj) {}
		PropertyMapping1(): _target(NULL) {}
	private:
		SrcT* _target;
	};
}
#endif 
