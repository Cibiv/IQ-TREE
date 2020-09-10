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
