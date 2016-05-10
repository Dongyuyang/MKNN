#include <SpatialIndex.h>
#include "include/common.h"

#define Entries_NUM 100

int main()
{
    int dimension = 2;

    IStorageManager *memfile =
        StorageManager::createNewMemoryStorageManager();

    StorageManager::IBuffer *file =
        StorageManager::createNewRandomEvictionsBuffer(*memfile, 10, false);

    id_type indexIdentifier;

    ISpatialIndex* tree =
        RTree::createNewRTree(*file, 0.7, Entries_NUM, Entries_NUM,
                                   dimension, SpatialIndex::RTree::RV_RSTAR,
                                   indexIdentifier);

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
    auto q = P[10];
    auto result = knn(tree, q, 5);
    cost.catch_time();

    cout << "N: " << P.size() << endl;
    cout << "Q:";
    displayCoordinates(q);cout << endl;
    //std::cout << "id: " << result << std::endl;
    displayPset(result);
    cout << "cpu cost is " << cost.get_cost(2) << " millisecond(s)" << endl;

}
