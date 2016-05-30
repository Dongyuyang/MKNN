#pragma once
#include "common.h"
#include "comcls.h"

class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
 public:
    MyQueryStrategy(const Point q, uint32_t k)
    {}

    void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext)
    {
        const INode *node = dynamic_cast<const INode*>(&entry);

        //internal node
        if(node->isIndex()){
            indexIO++;
            // for each child
            for(uint32_t n_child = 0; n_child < node->getChildrenCount(); n_child++){
                m_queue.push(new NNEntry(node->getChildIdentifier(n_child),0));
            }
        }
        //leaf node
        else{
            leafIO++;
            for(uint32_t n_child = 0; n_child < node->getChildrenCount(); n_child++){

            }
        }

        //priority_queue
        if(!m_queue.empty()){
            NNEntry* top = m_queue.top();
            m_queue.pop();
            nextEntry = top->m_id;
            hasNext = true;
        }else{
            hasNext = false;
        }

    }

 private:
    std::priority_queue<NNEntry*, vector<NNEntry*>, NNEntry::ascending> m_queue;


 public:
    uint32_t indexIO;
    uint32_t leafIO;
};
