#include "include/mbm.h"

#define CAPACITY 100 // index and leaf capacity
//#define FUN 1 // sum function
#define FUN 2 // max function
#define LOOPNUM 100

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		cerr << "Usage: " << argv[0] << " dim n area." << endl;
		return -1;
	}

	int dim = atol(argv[1]);
	int n = atol(argv[2]);
	double area = atof(argv[3]);
	if(dim <= 0)
	  {
		cerr << "Dimension should be larger than 0." << endl;
		return -1;
	}
	if(n <= 0)
	{
		cerr << "The number of query points should be larger than 0." << endl;
		return -1;
	}
	if(area <= 0 || area > 1)
	{
		cerr << "the area of query points should be in (0, 1]." << endl;
		return -1;
	}

	/*read static data set*/
	vector <Point> P;
	ifstream in("../data.ini");
	if(!in)
	{
		cerr << "Cannot open file data.ini.\n";
		return -1;
	}
 	P = readPoints(in, dim);
	uint32_t N = P.size();

	try {
		IStorageManager* memfile = StorageManager::createNewMemoryStorageManager();
		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*memfile, 10, false);
		id_type indexIdentifier;
		ISpatialIndex* tree = RTree::createNewRTree(*file, 0.7, CAPACITY, CAPACITY, dim, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
		id_type id = 0;
		for(uint32_t i = 0; i < N; ++i)
		{
			std::ostringstream os;
			os << P[i];
			std::string data = os.str();
			tree->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), P[i], id);
			id++;
		}

	for(uint32_t loop = 1; loop <= LOOPNUM; ++loop)
	{
		cout << "/**************** BEGIN " << loop << " ***************/" << endl;

		/*generate query set*/
		vector <Point> Q;
		//Q = genPoints(dim, n, area, loop);
		stringstream ss;
		ss << "../query/n" << n << "M" << area << "/loop" << loop;
		cout << ss.str().c_str() << endl;
		ifstream qin(ss.str().c_str());
		if(!qin)
		{
			cerr << "Cannot open query file";
			return -1;
		}
		Q = readPoints(qin, dim);
		
		/*************** BEGIN BF MBM method ******************/
		/* MBM method for finding ANN of Q */
		CATCH mbmcost;
		mbmcost.catch_time();

		Region M = getMBR(Q, dim, n);
		MyQueryStrategy qs(M, Q, n, FUN);
		tree->queryStrategy(qs);

		mbmcost.catch_time();
		cout << "MBM: cpu cost is " << mbmcost.get_cost(2) << " millisecond(s)" << endl;
		cout << "MBM: best_dist is " << qs.best_dist << endl;
		cout << "MBM: best_NN is ";
		displayCoordinates(qs.best_NN);
		cout << "MBM: leafIO = " << qs.mbm_leafIO << "; indexIO = " << qs.mbm_indexIO << endl << endl;
		/*************** END BF MBM method *******************/
		cout << "/**************** END " << loop << " ****************/" << endl << endl;

	} // end loop
		
		delete tree;
		delete file;
		delete memfile;
	}
	catch(Tools::Exception& e)
	{
		cerr << "*********ERROR**********" << endl;
		std::string s = e.what();
		cerr << s << endl;
		return -1;
	}
	catch(...)
	{
		cerr << "**********ERROR********" << endl;
		return -1;
	}

	return 1;
}

