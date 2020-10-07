//
//  terracetree.hpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#ifndef terracetree_hpp
#define terracetree_hpp

#include <stdio.h>
#include "tree/mtree.h"
#include "terracenode.hpp"

class TerraceTree: public MTree {

public:
    /**
     constructor
     */
    TerraceTree();
    
    /**
     constructor
     */
    //TerraceTree(SuperAlignment *alignment);
    
    /**
     destructor
     */
    ~TerraceTree();
    
    /**
     read the tree from the input file in newick format
     @param infile the input file file.
     @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(const char *infile, bool &is_rooted);
    
    /**
     read the tree from the ifstream in newick format
     @param in the input stream.
     @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(istream &in, bool &is_rooted);
    
    /**
     allocate a new node. Override this if you have an inherited Node class.
     @param node_id node ID
     @param node_name node name
     @return a new node
     */
    virtual Node* newNode(int node_id = -1, const char* node_name = NULL);
    
    /**
     allocate a new node. Override this if you have an inherited Node class.
     @param node_id node ID
     @param node_name node name issued by an interger
     @return a new node
     */
    virtual Node* newNode(int node_id, int node_name);
    
    /**
     copy the tree given a list of taxon names that should remain on the tree (not yet a 0-1 vector)
     */
    void copyTree_byTaxonNames(MTree *tree, vector<string> taxon_names);
    
    /**
     *  Clean all info about link neighbours and taxa
     */
    
    void cleanAllLinkINFO(TerraceNode *node = nullptr, TerraceNode *dad = nullptr);
    
};



#endif /* terracetree_hpp */
