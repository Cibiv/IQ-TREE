//
//  terracetree.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "terracetree.hpp"
#include "terracenode.hpp"

TerraceTree::TerraceTree(): MTree(){};
TerraceTree::~TerraceTree(){};


void TerraceTree::readTree(const char *infile, bool &is_rooted){
    MTree::readTree(infile, is_rooted);
};

void TerraceTree::readTree(istream &in, bool &is_rooted) {
    MTree::readTree(in, is_rooted);
}

Node* TerraceTree::newNode(int node_id, const char* node_name) {
    return (Node*) (new TerraceNode(node_id, node_name));
}

Node* TerraceTree::newNode(int node_id, int node_name) {
    return (Node*) (new TerraceNode(node_id, node_name));
}
