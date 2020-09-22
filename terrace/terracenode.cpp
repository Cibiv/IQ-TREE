//
//  terracenode.cpp
//  iqtree
//
//  Created by Olga on 09.09.20.
//

#include "terracenode.hpp"

TerraceNode::TerraceNode()
: Node()
{
    init();
}


TerraceNode::TerraceNode(int aid) : Node(aid)
{
    init();
}

TerraceNode::TerraceNode(int aid, int aname) : Node (aid, aname) {
    init();
}


TerraceNode::TerraceNode(int aid, const char *aname) : Node(aid, aname) {
    init();
}

void TerraceNode::init() {

}


void TerraceNode::addNeighbor(Node *node, double length, int id) {
    neighbors.push_back(new TerraceNeighbor(node, length, id));
}

TerraceNode::~TerraceNode()
{
}


void TerraceNeighbor::printInfo(Node *dad){
    
    cout<<"TerraceNeighbour information for node "<<this->node->id<<", a neighbour of node "<<dad->id<<":"<<endl;
    
    NeighborVec::iterator it1;
    vector<Node*>::iterator it2;
    
    int link_neighbors_size = link_neighbors.size();
    int taxa_to_insert_size = taxa_to_insert.size();
    
    int i=0;
    cout<<"There are "<<link_neighbors_size<<" link neighbors."<<endl;
    if(link_neighbors_size>0){
        for(it1=link_neighbors.begin(); it1<link_neighbors.end(); it1++){
            i++;
            cout<<"link_neighbor["<<i<<"]="<<(*it1)->node->id<<endl;
        }
        cout<<endl;
    }
    i=0;
    cout<<"There are "<<taxa_to_insert_size<<" taxa to insert."<<endl;
    if(taxa_to_insert_size>0){
        for(it2=taxa_to_insert.begin(); it2<taxa_to_insert.end(); it2++){
            i++;
            cout<<"taxa_to_insert["<<i<<"]="<<(*it2)->name<<endl;
        }
        cout<<endl;
    }
    
}
