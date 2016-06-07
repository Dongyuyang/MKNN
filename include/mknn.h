#pragma once
#include "common.h"
#include "comcls.h"

class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
 public:
    MyQueryStrategy(const Point query, uint32_t k_num)
        :q(query), k(k_num),best_dist(numeric_limits<double>::max()),
        indexIO(0),leafIO(0), counter(0)
    {}

    void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext)
    {
        const INode *node = dynamic_cast<const INode*>(&entry);

        //internal node
        if(node->isIndex()){
            indexIO++;
            // for each child
            for(uint32_t n_child = 0; n_child < node->getChildrenCount(); n_child++){
                IShape *shape;
                node->getChildShape(n_child, &shape);
                const Region *mbr = dynamic_cast<const Region *>(shape);
                if(mbr != NULL){
                    double distance = mbr->getMinimumDistance(q);
                    if(distance < best_dist)
                        m_queue.push(new NNEntry(node->getChildIdentifier(n_child),distance));
                }
                delete shape;
            }
        }
        //leaf node
        else{
            leafIO++;
            for(uint32_t n_child = 0; n_child < node->getChildrenCount(); n_child++){
                IShape *shape;
                node->getChildShape(n_child, &shape);
                const Region *mbr = dynamic_cast<const Region *>(shape);
                if(mbr != NULL){
                    Point p = Point(mbr->m_pLow, mbr->m_dimension);
                    double distance = q.getMinimumDistance(p);
                    if(distance < best_dist){
                        best_NN = p;
                        best_dist = distance;
                    }
                }
                delete shape;
            }
        }

        //priority_queue
        if(!m_queue.empty()){
            NNEntry* top = m_queue.top();
            m_queue.pop();
            nextEntry = top->m_id;
            hasNext = true;
            double dist = top->m_minDist;
            if(dist >= best_dist){
                counter++;
                best_KNN.push_back(best_NN);
                //best_dist = numeric_limits<double>::max();
                best_dist = top->m_minDist;
                if(counter == k)
                    hasNext = false;
            }
        }else{
            hasNext = false;
        }

    }

 private:
    std::priority_queue<NNEntry*, vector<NNEntry*>, NNEntry::ascending> m_queue;


 public:
    uint32_t indexIO;
    uint32_t leafIO;
    Point q;
    uint32_t k;
    double best_dist;
    Point best_NN;
    vector<Point> best_KNN;
    id_type best_NN_id;
    uint32_t counter;
    vector<double> candidates;

};
