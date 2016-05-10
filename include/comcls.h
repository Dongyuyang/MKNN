#pragma once

class NNEntry
{
	public:
		id_type m_id;
		double m_minDist;

		NNEntry(id_type id, double amindist) : m_id(id), m_minDist(amindist) {}
		~NNEntry() {}
		struct ascending : public std::binary_function<NNEntry*, NNEntry*, bool>
  	{
			bool operator()(const NNEntry* __x, const NNEntry* __y) const {return __x->m_minDist > __y->m_minDist; }
	  };
};


