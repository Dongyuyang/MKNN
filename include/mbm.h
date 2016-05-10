#pragma once

#include "common.h"
#include "comcls.h"

// the main idea of MBM method  (best first)
class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
	private:
		// list sorted on amindist<N_j, Q>
		std::priority_queue<NNEntry*, std::vector<NNEntry*>, NNEntry::ascending> m_queue; 
		Region m_M;
		vector <Point> m_Q;
		uint32_t m_n;
		uint32_t m_f; // aggregate function

	public:
		Point best_NN;
		double best_dist;
		uint32_t mbm_indexIO;
		uint32_t mbm_leafIO;

	public:
		MyQueryStrategy(const Region &M, const vector <Point> &Q, const uint32_t n, const uint32_t f)
			: m_M(M), m_Q(Q), m_n(n), m_f(f), 
			  best_dist(numeric_limits<double>::max()),
        mbm_indexIO(0), mbm_leafIO(0)	{} 

		void getNextEntry(const IEntry &entry, id_type& nextEntry, bool& hasNext)
		{ 
			const INode* node = dynamic_cast<const INode*>(&entry);
			// If Node is an intermediate node
			if(node->isIndex())
			{
				mbm_indexIO++;

				// for each entry N_j of Node
				for(uint32_t cChild = 0; cChild < node->getChildrenCount(); cChild++)
				{
					IShape* pS;
					node->getChildShape(cChild, &pS);
					const Region *pr = dynamic_cast<const Region *>(pS);
					if(pr != NULL)
					{
			  		//double amindist3 = n * pr.getMinimumDisntace(M);
						double amindist3 = getAmindist(*pr, m_M, m_n, m_f); // amindist(N_j, M
						if(amindist3 < best_dist) // amindist(N_j, M) < best_dist   
						{
							double amindist2 = getAmindist(*pr, m_Q, m_n, m_f);
							//cout << "amindist3 = " << amindist3 << "; amindist2 = " << amindist2 << "; best_dist = " << best_dist << endl;
							if(amindist2 < best_dist) 	// amindist(N_j, Q) < best_dist
							{
								// insert N_j in a list sorted on amindist(N_J, Q);
								m_queue.push(new NNEntry(node->getChildIdentifier(cChild), amindist2)); 
		  				}
							else { /* amindist(N_j, Q) >= best_dist */ }
						} 
						else { /* amindist(N_j, M) >= best_dist */ }
					}
					delete pS;
				}
				//cout << endl; 
		  }
			else // if Node is a leaf node
			{
				mbm_leafIO++; 

				// for each data point pj in Node 
				for(uint32_t cChild = 0; cChild < node->getChildrenCount(); cChild++)
				{
					IShape* pS;
					node->getChildShape(cChild, &pS);
					const Region *pr = dynamic_cast<const Region *>(pS);
					if(pr != NULL) // Region
					{  
						uint32_t dim = pr->m_dimension;
						double * pCoords = pr->m_pLow;
						Point p = Point(pCoords, dim);
						double amindist3 = getAmindist(p, m_M, m_n, m_f); // amindist(pj, M)]
						if (amindist3 < best_dist) // amindist(pj, M) < best_dist
						{
							double adist = getAdist(p, m_Q, m_n, m_f); 
							if(adist < best_dist) // adist(pj, Q) < best_dist
							{	
								best_NN = p;						
								best_dist = adist;
							}
							else { /* adist(pj, Q) >= best_dist */ }
						} 
						else { /* amindist(pj, M) >=  best_dist */ }
					} 
					else { /*	Region pr is NULL */ }
					delete pS;
				} 
			} // end else if Node is a leaf node

		//	while(! m_queue.empty())
		  if (! m_queue.empty())
			{
				// get next entry N_j from list]
				NNEntry* pFirst = m_queue.top();
				m_queue.pop();
				nextEntry =  pFirst->m_id;
				hasNext = true;
				double amindist2 = pFirst->m_minDist; 
				if (amindist2 < best_dist) // amindist(N_j, Q) < best_dist 
				{
					// Recursion
					//break;
				}
				else
				{ /* amindist(N_j, Q) >= best_dist */
					hasNext = false;
					//break;
				}
			}	
			//if(m_queue.empty())
			else 
				hasNext = false;
		}
};
