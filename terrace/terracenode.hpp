//
//  terracenode.hpp
//  iqtree
//
//  Created by Olga on 09.09.20.
//

#ifndef terracenode_hpp
#define terracenode_hpp

#include <stdio.h>
#include "node.h"

/**
 A neighbor in a tree, which belongs to some terrace
 
 @author Olga Chernomor <olga.chernomor@univie.ac.at>
 */
class TerraceNeighbor : public Neighbor {
    
    friend class TerraceNode;
    friend class TerraceTree;
    
public:
    /**
     construct class with a node and length
     @param anode the other end of the branch
     @param alength length of branch
     */
    TerraceNeighbor(Node *anode, double alength) : Neighbor(anode, alength) {
    }
    
    /**
     construct class with a node and length
     @param anode the other end of the branch
     @param alength length of branch
     @param aid branch ID
     */
    TerraceNeighbor(Node *anode, double alength, int aid) : Neighbor(anode, alength, aid) {
    }
    
    // TODECIDE: for simplicity you might want to have just vector of nodes, oder?
    
    /**
     vector of size m (m = #partitions) for representative tree and 1 for a back map from induced partition tree to the representative (correct?)
     */
    //NeighborVec link_neighbors;
    //NeighborVec link_neighbors_taxa;
    
};

/**
 Node of a tree on a terrace
 
 @author Olga Chernomor <olga.chernomor@univie.ac.at>
 */
class TerraceNode : public Node
{
    friend class TerraceTree;
    
public:
    /**
     constructor
     */
    TerraceNode();
    
    /**
     constructor
     @param aid id of this node
     */
    TerraceNode(int aid);
    
    /**
     constructor
     @param aid id of this node
     @param aname name of this node
     */
    TerraceNode(int aid, int aname);
    
    /**
     constructor
     @param aid id of this node
     @param aname name of this node
     */
    TerraceNode(int aid, const char *aname);
    
    /**
     initialization
     */
    void init();
    
    /**
     add a neighbor
     @param node the neighbor node
     @param length branch length
     @param id branch ID
     */
    virtual void addNeighbor(Node *node, double length, int id = -1);
    
    ~TerraceNode();
    
};


#endif /* terracenode_hpp */
