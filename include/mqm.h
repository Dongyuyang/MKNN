#pragma once

#include "common.h"

#define NEXTKNN 100
#define ONEIO 1000  // 1 I/O can read or wirte 1000 of uint32_t id

uint32_t mqm_indexIO = 0;
uint32_t mqm_leafIO = 0;

/*
 * Finding kNN for Point q in Points Set P
 */
vector <Point> retriveKNNs(ISpatialIndex* tree, const Point & q, uint32_t k) {
	vector <Point> kNNs;
	try {
		MyVisitor vis;
		tree->nearestNeighborQuery(k, q, vis); // find kNN for Point q
		kNNs = vis.m_kNNs;
		
		// leaf and index IO for finding k nearest neighbor
		mqm_indexIO += vis.m_indexIO;  
		mqm_leafIO += vis.m_leafIO;

		return kNNs;
	}
	catch (Tools::Exception& e) {
		cerr << "******ERROR******" << endl;
		std::string s = e.what();
		cerr << s << endl;
		return kNNs;
	}
	catch (...) {
		cerr << "******ERROR******" << endl;
		cerr << "other exception" << endl;
		return kNNs;
	}
}

/*
 * Update global threshold T
 */
double updateT(double t[], uint32_t n, uint32_t f)
{
	double newT = 0.0;
	switch(f)
	{
		case 1:  // f = sum
			{
				for (uint32_t i = 0; i < n; ++i)
					newT += t[i];
				break;
			}
		case 2:  // f = max
			{
				for (uint32_t i = 0; i < n; ++i)
					newT = max(newT, t[i]);
				break;
			}
		case 3:  // f = min
			{
				newT = numeric_limits<double>::max();
				for (uint32_t i = 0; i < n; ++i)
					newT = min(newT, t[i]);
				break;
			}
	}
	return newT;
}

/************************************************************************************/
void MQM(ISpatialIndex* tree, const vector <Point> & Q, const uint32_t n, const uint32_t f, double epsilon)
{
	CATCH mqmCatch;
	CATCH knnCatch;
	double knnCost = 0.0;
	double all_knnCost = 0.0;
	mqmCatch.catch_time();
	knnCatch.catch_time();

	uint32_t k = NEXTKNN;
	vector <Point> QkNNs[n];
	for(uint32_t i = 0; i < n; ++i)
	{  // find kNNs (k = NEXTKNN)
	  QkNNs[i] = (retriveKNNs(tree, Q[i], k));
	}

	knnCatch.catch_time();
	all_knnCost = knnCost = knnCatch.get_cost(2);
 
	/* T: global threshold; best_dist: aggregate distance of the current NN */
	/* Initialization */
	//double t[n];
	//for(uint32_t i = 0; i < n; i++) t[i] = 0;
	double T = 0.0, best_dist = numeric_limits<double>::max(); 
	Point best_NN;  	
	uint32_t times = 0;
	while (best_dist > (1+epsilon) * T) 
	{
		if(times >= k)
	 	{
		  knnCatch.catch_time();

			mqm_indexIO = 0; mqm_leafIO = 0; // IO cost
			k += NEXTKNN;
			for(uint32_t i = 0; i < n; ++i)
			{ // find kNNs (k += NEXTKNN)
				QkNNs[i] = retriveKNNs(tree, Q[i], k);
			}

			knnCatch.catch_time();
			knnCost = knnCatch.get_cost(2);
			all_knnCost += knnCost;
		} 

		for(uint32_t i = 0; i < n; ++i){  // select the next query point qi
			Point pj = QkNNs[i][times]; // get the next nearest neighbor pj of qi;
			
			//f= sum
			double newti = pj.getMinimumDistance(Q[i]);  // ti = |pj qi|
			//T = T - t[i] + newti; // update T, T=f(t1, t2, ..., tn)
			//t[i] = newti; 
			
			//f = max
			T = max(T, newti);

			double tmp_adist = getAdist(pj, Q, n, f); // compute adist(pj, Q)
			// Update current ANN of Q
			if(tmp_adist < best_dist){
				best_NN = pj; best_dist = tmp_adist;
			} // Update current ANN of Q
		
			if(best_dist <= (1+epsilon) * T) break;  
		}   
		times++;
	}

	mqmCatch.catch_time();
	// show the result information
	cout << epsilon << "-MQM: maxKNN for MQM is:" << times << endl;
	cout << epsilon << "-MQM: cpu cost is " << mqmCatch.get_cost(2) - all_knnCost + knnCost << " millisecond(s)" << endl;
	cout << epsilon << "-MQM: best_dist is " << best_dist << endl;
	cout << epsilon << "-MQM: best_NN is ";
  displayCoordinates(best_NN);
	cout << epsilon << "-MQM: leafIO = " << mqm_leafIO << "; indexIO = " << mqm_indexIO << endl << endl; // IO cost
}

/*
 * using the main idea of MQM to find ANN of Q with dist-matrix index, (Q in P)
 * @param: const vector <Point> Q: query set
 * @param: const uint32_t n: the number of Q
 * @param: const vector <Point> P: static data set
 * @param: const uint32_t N: the number of P
 * @param: uint32_t** dmindex:  dist-matrix index of Q
 * @param: uint32_t maxK:   the length of dist-matrix index
 * @param: const uint32_t f:  aggregate function (1: sum; 2: max; 3: min)
 * @param: const double minadist: lower bound
 * @param: const double err: eplison (error rate)
 */
Point ADM(const vector <Point> & Q, const uint32_t n,
		const vector <Point> & P, const uint32_t N,
		uint32_t** dmindex, const uint32_t maxK, const uint32_t f,
		const double minadist, const double err)
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

	CATCH admCatch;
	admCatch.catch_time();
 
	/* T: global threshold; best_dist: aggregate distance of the current NN */
	/* Initialization */
  //double t[n];
  //for(uint32_t i = 0; i < n; i++) t[i] = 0;
	double T = 0.0, best_dist = numeric_limits<double>::max(); 
	Point best_NN;  	
	uint32_t times = 0;
	while (T < best_dist) 
	{
		//if (best_dist <= (1+ err)* minadist) break;

		for(uint32_t i = 0; i < n; ++i){  // select the next query point qi
			uint32_t id = dmindex[i][times]; // get the next nearest neighbor pj of qi;
			//f = sum
			double newti = P[id].getMinimumDistance(Q[i]);  // ti = |pj qi|
			//T = T - t[i] + newti; // update T, T=f(t1, t2, ..., tn)
			//t[i] = newti; 

			//f = max
			T = max(T, newti);

			if(!visit[id])  // p_ix hasn't been visited
			{
				double tmp_adist = getAdist(P[id], Q, n, f); // compute adist(pj, Q)
				// Update current ANN of Q
				if(tmp_adist < best_dist){
					best_NN = P[id]; best_dist = tmp_adist;
				} // Update current ANN of Q
	
				if(T >= best_dist) break;  

			  visit[id] = true;  // p_ix is visited
			} 
			else {}
		} // end for q'_i
		times++;
	}

	admCatch.catch_time();
	cout << "ADM: cpu cost is " << admCatch.get_cost(2) << " millisecond(s)" << endl;
	cout << "ADM: best_NN for Q' is ";
	displayCoordinates(best_NN);
	cout << "ADM: kNN for read dist-matrix is: times * n = " << times * n << endl; 

	return best_NN;
}
