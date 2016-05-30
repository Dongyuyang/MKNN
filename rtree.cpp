#include "include/mknn.h"

#define Entries_NUM 100

int main(int argc, char* argv[])
{
    /*default params*/
    int dimension = 2;

    /*input params*/
    int k = atoi(argv[1]);

    /*Init Rtree*/
    IStorageManager *memfile = StorageManager::createNewMemoryStorageManager();
    StorageManager::IBuffer *file = StorageManager::createNewRandomEvictionsBuffer(*memfile, 10, false);
    id_type indexIdentifier;
    ISpatialIndex* tree2 = RTree::createNewRTree(*file, 0.7, Entries_NUM, Entries_NUM,
                                                 dimension, SpatialIndex::RTree::RV_RSTAR,
                                                 indexIdentifier);
    auto tree = RTree::loadRTree(*file,indexIdentifier);

    /*read data insert rtree*/
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

    /*run*/
    CATCH cost;
    cost.catch_time();
    auto q = P[20];
    MyQueryStrategy qs(q,5);
    tree->queryStrategy(qs);
    cout << "MY: leafIO = " << qs.leafIO << "; indexIO = " << qs.leafIO << endl;
    cost.catch_time();

    cout << "my nn id: " << qs.best_NN_id << endl;
    displayCoordinates(qs.best_NN);
    auto r_id = nearestNeighbor(tree, q);
    cout << "r_id: " << r_id << endl;

    //CATCH costb;
    // costb.catch_time();
    // auto result = mrnn(tree,q,k);
    // costb.catch_time();
    //auto re = knn(tree,q,k);

    cout << "N: " << P.size() << endl;
    cout << "Q:"; displayCoordinates(q); cout << endl;
    cout << "cpu cost is " << cost.get_cost(2) << " millisecond(s)" << endl;
    //cout << "mrnn cost is " << costb.get_cost(2) << " millisecond(s)" << endl;

    /*release*/
    delete tree;
    delete file;
    delete memfile;

}
