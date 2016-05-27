#pragma once
#pragma warning(disable:4786)

#include <cstring>
#include <stdint.h>
#include <cmath>
#include <limits>
//#include <SpatialIndex.h>
#include "spatialindex/SpatialIndex.h"
#include "visitor.h"       // for using rtree to find nearest neighbor
#include "catch.h"         // catch runing time

using namespace SpatialIndex;
using namespace std;


void displayData(const double * pData, const uint32_t dimension)
{
	if (pData == NULL || dimension <= 0)
	{
		cout << "Error: data is NULL or dimension is lower than 1.\n";
	}
	else 
	{
		for(uint32_t i = 0; i < dimension; ++i)
		{
			if (i == dimension - 1)
				cout << pData[i] << endl;
			else
				cout << pData[i] << " ";
		}
	}
}

/**************** Begin: Operatioins of Data Point **************************/
/*
 * display the coordinates of data point p
 * @param: Point vp
 * @return: void
 */
void displayCoordinates(const Point &vp)
{
	for(uint32_t i = 0; i < vp.m_dimension; ++i)
	{
		if( i == vp.m_dimension -1 )
			cout << vp.m_pCoords[i] << endl;
		else
			cout << vp.m_pCoords[i] << " ";
	}
}

/*
 * show the cordinates of data points
 * @param: vector <Point> pset
 * @return: void
 */
void displayPset(const vector <Point>& pset)
{
	uint32_t pcount = pset.size();
	for(uint32_t i = 0; i < pcount; i++)
	{
		displayCoordinates(pset[i]);
	}
}

/*
 * generate data set into [0, 1]
 * @param: const uint32_t dimension
 * @param: const uint32_t count
 * @param: const double area, the ratio of the area of data point to all area
 * @return: vector <Point> pset
 */
vector <Point> genPoints(const uint32_t dimension, 
		const uint32_t count, const double area, uint32_t seed)
{	
	double power = (double) 1 / dimension;
	double gap = pow(area, power);

	vector<Point> pset;
	Tools::Random rnd = Tools::Random(seed, 0xD31A);
	double measure[dimension];  // lower bound
	for(uint32_t j = 0; j < dimension; j++)
	{
		if(gap == 1){
			measure[j] = 0.0;
			continue;
		}

		measure[j] = rnd.nextUniformDouble();
		while(measure[j] > (1 - gap))
			measure[j] = rnd.nextUniformDouble();
	}

	for(uint32_t i = 0; i < count; i++)
	{
		double * coords = new double[dimension];
		for(uint32_t j = 0; j < dimension; j++)
		{
			coords[j] = rnd.nextUniformDouble();
			while (coords[j] < measure[j] || (coords[j] - measure[j]) > gap)
				coords[j] = rnd.nextUniformDouble();
		}
		
		Point p = Point(coords, dimension);
		pset.push_back(p);
		//displayCoordinates(p);
	} 

	return pset;
}

/*
 * wirte data set into file
 */
bool writePoints(const vector <Point> &pset, const char* filename)
{
	ofstream out(filename);
	if(!out)
	{
		cerr << "Cannot open file.\n";
		return false;
	}

	uint32_t count = pset.size();
	uint32_t dimension = pset[0].m_dimension;

	for(uint32_t i = 0; i < count; ++i)
	{
		for(uint32_t j = 0; j < dimension; ++j)
		{
			out << pset[i].m_pCoords[j] << " ";
		}
		out << endl;
	}
	out.close();
	return true;
}

/*
 * read data set from file
 * @param: ifstream& fin
 * @param: const uint32_t dimension
 * @return: vector <Point> pset
 */
vector <Point> readPoints(ifstream& fin, const uint32_t dimension)
{
	vector <Point> pset;

	double pdata[dimension];
    uint32_t id = 1;
	while(fin)
	{
		for(uint32_t i = 0; i < dimension; i++)
			fin >> pdata[i];
		if(! fin.good()) continue;
		Point vp(pdata, dimension);
        vp.m_id = id;
		pset.push_back(vp);
        id++;
	}
	return pset;
}

/*
 * find centroid point of point set
 * @return: Point centroid
 */
Point findCentroid(const vector <Point> &pset, 
		const uint32_t dimension, const uint32_t count)
{
	double coords[dimension];
	for(uint32_t j = 0; j < dimension; j++)
		coords[j] = 0.0;

	for(uint32_t j = 0; j < dimension; j++)
	{
		for(uint32_t i = 0; i < count; i++)
			coords[j] += pset[i].m_pCoords[j];
	}

	for(uint32_t j = 0; j < dimension; j++)
		coords[j] = coords[j] / count;

	Point centroid(coords, dimension);
	return centroid;
}

/*
 * find centroid point of point set
 * @return: Point centroid
 */
Point findCentroid(const vector <uint32_t> &nnId, const vector <Point>& P,
		const uint32_t dimension, const uint32_t count)
{
	vector <Point> pset;
	for(uint32_t i = 0; i < count; ++i)
		pset.push_back(P[nnId[i]]);
	Point centroid = findCentroid(pset, dimension, count);
	return centroid;
}

/*
 * find lower bound of MBR
 */
Point findLowerBound(const vector <Point> &pset, 
		const uint32_t dimension, const uint32_t count)
{
	double coords[dimension];
 	for(uint32_t j = 0; j < dimension; j++)
		coords[j] = numeric_limits<double>::max();

	for(uint32_t i = 0; i < count; i++)
	{
		for (uint32_t j = 0; j < dimension; j++)
			coords[j] = min(coords[j], pset[i].m_pCoords[j]);
	}

	Point lowerBound(coords, dimension);
	return lowerBound;
}

/*
 * find high bound of MBR
 */
Point findHighBound(const vector <Point> &pset, 
		const uint32_t dimension, const uint32_t count)
{
	double coords[dimension];
	for(uint32_t j = 0; j < dimension; j++)
		coords[j] = 0.0;

	for(uint32_t i = 0; i < count; i++)
	{
		for(uint32_t j = 0; j < dimension; j++)
			coords[j] = max(coords[j], pset[i].m_pCoords[j]);
	}

	Point highBound(coords, dimension);
	return highBound;
}

Region getMBR(const vector <Point> &pset,
		const uint32_t dimension, const uint32_t count)
{
	Point lowbound = findLowerBound(pset, dimension, count);
	Point highbound = findHighBound(pset, dimension, count);
	Region MBR = Region(lowbound, highbound);
	return MBR;
}

/**************** End: Operatioins of Data Point **************************/

/**************** Begin: Operatioins of Dist-Matrix **************************/
/*
 * read dist-matrix of point p from file 
 * @param: const uint32_t id, id of point p, 
 *         and the corresponding filename is id
 * @param: const uint32_t maxK
 * @return: uint32_t * index
 */
uint32_t* readDMIndex(const uint32_t id, const uint32_t maxK)
{
	stringstream ss;
	ss << id;
	const char * filename = ss.str().c_str();
	ifstream in(filename, ios::in | ios::binary);
	if(!in)
	{
		cout << "Cannot open file.\n";
		return NULL;
	}

	uint32_t * index = new uint32_t[maxK];
	in.read((char *) index, (sizeof index) * maxK);
	in.close();

	return index;
}

/*
 * display dist-matrix of point p
 * @param: const uint32_t [] index
 * @param: const uint32_t masK
 * @return: void
 */
void displayDMIndex(const uint32_t * index, const uint32_t maxK)
{
	for(uint32_t i = 0; i < maxK; ++i)
		cout << index[i] << " ";
	cout << endl;
}

/**************** End: Operatioins of Dist-Matrix **************************/
/*
 * using R-tree to find nearest neighbor of query point q
 * @return: id of nearest neighbor point
 */
uint32_t nearestNeighbor(ISpatialIndex* tree, const Point &q)
{
	nnVisitor nnvis;
	uint32_t firstNN = 1;
	tree->nearestNeighborQuery(firstNN, q, nnvis); // find nearest neighbor for q

	// index and leaf IO for fining nearest neighbor
	cout << "indexIO: "  << nnvis.m_indexIO << "; leafIO: " << nnvis.m_leafIO << endl;

	return nnvis.getIdentifierId();
}

vector<Point> knn(ISpatialIndex* tree, const Point &q, uint32_t k)
{
    MyVisitor knnvis;
    tree->nearestNeighborQuery(k,q,knnvis);
    cout << "indexIO: " << knnvis.m_indexIO << "; leafIO: " << knnvis.m_leafIO << endl;
    return knnvis.m_kNNs;
}

vector<Point> mrnn(ISpatialIndex* tree, const Point &q, uint32_t k)
{
    vector<Point> result;
    auto kn = knn(tree,q,k);
    for(auto nn : kn){
        auto t_knn = knn(tree,nn,k);
        int counter = 0;
        for(auto tk : t_knn){
            if(tk == q){
                result.push_back(nn);
                /*counter++;
                if(counter == 2){
                    result.push_back(nn);
                    }*/
            }
        }
    }
    return result;
}

/*
 * using R-tree to find nearest neighbor for every query point q in qset
 * @return: return ids of the corresponding nearest neighbor points
 */
vector <uint32_t> nearestNeighborSet(ISpatialIndex* tree, const vector <Point> &qset, uint32_t n)
{
	uint32_t indexIO = 0, leafIO = 0;

	//IDs of nearest neighbor point for all query point
	vector <uint32_t> nnIDs;
	for(uint32_t i = 0; i < n; ++i)
	{
		//uint32_t id = nearestNeighbor(tree, qset[i]);
		nnVisitor nnvis;
		tree->nearestNeighborQuery(1, qset[i], nnvis); // find nearest neighbor for qi
		uint32_t id = nnvis.getIdentifierId();
		nnIDs.push_back(id);

    /* index and leaf IO for finding NN of all qi with R-tree*/
		indexIO += nnvis.m_indexIO; leafIO += nnvis.m_leafIO;
	}
	return nnIDs;
}

/*
 * calculate aggregate distance of p to Q with monotonic function f
 */
double getAdist(const Point & p, const vector <Point>& qset, const uint32_t n, const uint32_t f)
{
	double adist = 0.0;
	switch(f)
	{
		case 1: // f = sum
			{
				for(uint32_t i = 0; i < n; ++i)
					adist += p.getMinimumDistance(qset[i]);
				break;
			}
		case 2: // f = max
			{
				for(uint32_t i = 0; i < n; ++i)
					adist = max(adist, p.getMinimumDistance(qset[i]));
				break;
			}
		case 3: // f = min
			{
				adist = numeric_limits<double>::max();
				for(uint32_t i = 0; i < n; ++i)
					adist = min(adist, p.getMinimumDistance(qset[i]));
				break;
			}
	}
	return adist;
}

// calculate amindist(N_j, M)
double getAmindist(const Region& pr, const Region& M, const uint32_t n, uint32_t f)
{
	double amindist = pr.getMinimumDistance(M);
	if (f == 1) // sum
		return amindist * n;
	else  // f == 2 (max) or f == 3 (min)
		return amindist;
}

// calculate amindist(N_j, Q)
double getAmindist(const Region& pr, const vector <Point>& qset, const uint32_t n, uint32_t f)
{
	double amindist = 0.0;
	switch(f)
	{
		case 1:  // f = sum
			{
				for(uint32_t i = 0; i < n; ++i)
					amindist += pr.getMinimumDistance(qset[i]);
				break;
			}
		case 2:  // f = max
			{
				for(uint32_t i = 0; i < n; ++i)
					amindist = max(amindist, pr.getMinimumDistance(qset[i]));
				break;
			}
		case 3:  // f = min
			{
				amindist = numeric_limits<double>::max();
				for(uint32_t i = 0; i < n; ++i)
					amindist = min(amindist, pr.getMinimumDistance(qset[i]));
				break;
			}
	}
	return amindist;
}

// calculate amindist(pj, M)
double getAmindist(const Point& p, const Region& M, const uint32_t n, uint32_t f)
{
	double amindist = p.getMinimumDistance(M);
	if (f == 1) // sum
		return amindist * n;
	else // max or min
		return amindist;
}

/*********************************************************************************/
/*
 * @param: uint32_t vpID, id of nearest neighbor point vp for centroid point vq of Q
 * @return: uint32_t ID, id of epsilon ANN
 */
uint32_t epsilonANN(const vector <Point> &qset, const uint32_t n,
	 	const vector <Point> &pset, const uint32_t N,
		const uint32_t vpID, 
		uint32_t** dmindex, const uint32_t maxK, const uint32_t f)
{
	/*
	 * if pi has visited, visit[i] = true
	 * else visit[i] = false;
	 */
	vector <bool> visit;
	for(uint32_t i = 0; i < N; ++i)
	{
		visit.push_back(false);
	}

	CATCH cost;
	cost.catch_time();

	visit[vpID] = true;
	//best_dist: adist of vp to Q as upper bound
	double best_dist = getAdist(pset[vpID], qset, n, f);
	double temp_dist = 0.0;
	uint32_t best_NN_id = vpID;

	for(uint32_t i = 0; i < n; ++i)
	{
		for(uint32_t j = 0; j < maxK; ++j)
		{
			uint32_t ix = dmindex[i][j];
			if(ix == vpID) break; // the second loop terminates when visited vp

			if (!visit[ix]) // p_ix hasn't been visited
			{
				temp_dist = getAdist(pset[ix], qset, n, f);
				if(temp_dist < best_dist){
					best_dist = temp_dist;
					best_NN_id = ix; 
				}
				visit[ix] = true;  // p_ix is visited
			} 
			else {/* p_ix has been visited */ }
		} // end for q'_i
	} // end for Q

	cost.catch_time();
	cout << "approxiamte vp-ANN method: cpu cost is " << cost.get_cost(2) << " millisecond(s)" << endl;

	return best_NN_id;
}

/*
 * find ANN of Q' as epsilon ANN of Q
 */
uint32_t epsilonANN(const vector <uint32_t> &newQ_IDs, const uint32_t n, 
		const vector <Point>& pset, const uint32_t N,
		const uint32_t vpID,
		uint32_t** dmindex, const uint32_t maxK, const uint32_t f)
{
	vector <Point> Q_pi;
	for(uint32_t i = 0; i < n; ++i)
		Q_pi.push_back(pset[newQ_IDs[i]]);

	return epsilonANN(Q_pi, n, pset, N, vpID, dmindex, maxK, f);
}

/*
 * brute method for finding ANN of Q
 */
uint32_t brute_ANN(const vector <Point> &Q, const uint32_t n,
		const vector <Point> &P, const uint32_t N,
		const uint32_t f)
{
	uint32_t ANNid = 0;
	double best_dist = getAdist(P[ANNid], Q, n, f);
	double temp_dist = 0.0;
	for(uint32_t i = 1; i < N; ++i)
	{
		temp_dist = getAdist(P[i], Q, n, f);
		if(temp_dist < best_dist)
		{
			best_dist = temp_dist;
			ANNid = i;
		}
	}
	return ANNid;
}
