//#include "include/spatialindex/SpatialIndex.h"
#include "include/common.h"

#define Entries_NUM 100

int main(int argc, char* argv[])
{
    int dimension = 2;
    int k = atoi(argv[1]);

    IStorageManager *memfile =
        StorageManager::createNewMemoryStorageManager();

    StorageManager::IBuffer *file =
        StorageManager::createNewRandomEvictionsBuffer(*memfile, 10, false);

    id_type indexIdentifier;

    ISpatialIndex* tree2 =
        RTree::createNewRTree(*file, 0.7, Entries_NUM, Entries_NUM,
                                   dimension, SpatialIndex::RTree::RV_RSTAR,
                                   indexIdentifier);
    auto tree = RTree::loadRTree(*file,indexIdentifier);

    std::vector<Point> P;
    std::ifstream in("data/data2.ini");
    P = readPoints(in, dimension);
    id_type id = 0;
    for(uint32_t i = 0 ; i < P.size(); i++){
        std::ostringstream os;
        os << P[i];
        std::string data = os.str();
        tree->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), P[i], id);
        id++;
    }

    CATCH cost;
    cost.catch_time();
    //int result = nearestNeighbor(tree, P[1]);
    auto q = P[20];
    std::vector<double> t_q = {2,2};
    Point t_qq(&t_q[0],2);
    //auto result = knn(tree, q, 5);
    auto result = mrnn(tree,q,k);

    cost.catch_time();

    CATCH cost2;
    cost2.catch_time();
    auto t = knn(tree,q,k);
    cost2.catch_time();

    cout << "N: " << P.size() << endl;
    cout << "Q:";
    displayCoordinates(t_qq);
    cout << endl;
    //std::cout << "id: " << result << std::endl;
    displayPset(result);
    cout << "cpu cost is " << cost.get_cost(2) << " millisecond(s)" << endl;
    cout << "cpu2 cost is " << cost2.get_cost(2) << " millisecond(s)" << endl;


    delete tree;
    delete file;
    delete memfile;

}
