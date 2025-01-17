#pragma once

#include <SpatialIndex.h>

using namespace SpatialIndex;
using namespace std;

#define NUMITER 50

//Calculate the first candidate median as the geometric mean
Point candMedian(const vector <Point> &pset,
		const uint32_t dimension, const uint32_t N)
{
	double coords[dimension];
	for(uint32_t j = 0; j < dimension; ++j)
		coords[j] = 0.0;

	for(uint32_t i = 0; i < N; ++i)
	{
		for(uint32_t j = 0; j < dimension; ++j)
			coords[j] += pset[i].m_pCoords[j];
	}

	for(uint32_t j = 0; j < dimension; ++j)
		coords[j] = coords[j] / N;

	Point mean(coords, dimension);
	return mean;
}

/* Provides the denominator of the weiszfeld algorithm 
 * depending on whether you are adjusting the candiate x or y
 */
double numersum(const Point &p, const Point &q)
{
	return 1 / p.getMinimumDistance(q);
}

/* rovides the denominator of the weiszfeld algorithm */
double denomsum(const Point &p, const vector <Point> &pset, const uint32_t N)
{
	double temp = 0.0;
	for (uint32_t i = 0; i < N; ++i)
	{
		temp += 1 / p.getMinimumDistance(pset[i]);
	}
	return temp;
}

/* This function calculates the sum of linear distances from 
 * the current candidate median to all points in the data set, 
 * as such it is the objective function we are minimising.
 */
double objfunc(const Point &p, const vector <Point> &pset, const uint32_t N)
{
	double temp = 0.0;
	for(uint32_t i = 0; i < N; ++i)
		temp += p.getMinimumDistance(pset[i]);

	return temp;
}

/* minimise the objective function */
Point findMedian(const vector <Point> &pset, const uint32_t dimension, const uint32_t N)
{
	Point median = candMedian(pset, dimension, N);
	double tmpadist = denomsum(median, pset, N);
	double nextCoords[dimension];
	for(uint32_t i = 0; i < NUMITER; ++i)
	{
		double denom = denomsum(median, pset, N);
		for(uint32_t j = 0; j < dimension; ++j)
		{
			nextCoords[j] = 0.0; 
		}

		for(uint32_t i = 0; i < N; ++i)
		{
			for(uint32_t j = 0; j < dimension; ++j)
			{
				nextCoords[j] += (pset[i].m_pCoords[j] * numersum(median, pset[i])) / denom;
			}
		}

		for(uint32_t j = 0; j < dimension; ++j)
			median.m_pCoords[j] = nextCoords[j];

		double adist = objfunc(median, pset, N);
		if (tmpadist == adist) break;
		else tmpadist = adist;
		//cout << "i = " << i << "; adist = " << adist << endl;
	}
	return median;
}


/******************** for max function ***************************/
uint32_t findFurthestNeighbor(const vector <Point> &pset, const uint32_t N, const Point &p)
{
	double maxdist = 0;	uint32_t maxId = 0;
	for (uint32_t i = 0; i < N; i++)
	{
		double tmpdist = p.getMinimumDistance(pset[i]);
		if(maxdist < tmpdist)
		{
			maxdist = tmpdist;
			maxId = i;
		}
	}
	return maxId;
}

double cal_radius_MEB(const vector <Point> &pset, const uint32_t dimension, const uint32_t N, const double epsilon)
{
	uint32_t alpha = findFurthestNeighbor(pset, N, pset[0]);
	uint32_t beta = findFurthestNeighbor(pset, N, pset[alpha]);
	double pCoords[dimension];
	for(uint32_t i = 0; i < dimension; i++)
	{
		pCoords[i] = (pset[alpha].m_pCoords[i] + pset[beta].m_pCoords[i]) / 2;
	}
	Point center(pCoords, dimension);  // center of MEB
	double gamma = pow(pset[alpha].getMinimumDistance(pset[beta]) / 2, 2);
	uint32_t kappa = findFurthestNeighbor(pset, N, center);
	double delta = pow(center.getMinimumDistance(pset[kappa]), 2) / gamma - 1;

	while (delta > pow((1+epsilon), 2) - 1)
	{
		double lambda = delta / (2 * (1 + delta));
		for(uint32_t i = 0; i < dimension; i++)
		{
		  center.m_pCoords[i] = (1 - lambda) * center.m_pCoords[i] + lambda * pset[kappa].m_pCoords[i];
		}
		kappa = findFurthestNeighbor(pset, N, center);
		gamma = gamma * (1 + pow(delta, 2) / (4* (1 + delta)));
	  delta = pow(center.getMinimumDistance(pset[kappa]), 2) / gamma - 1;
	}

	return gamma;
}
