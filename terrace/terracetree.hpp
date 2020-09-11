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

class TerraceTree: public MTree, public vector<MTree*> {

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

    MTree *master_tree;
    vector<MTree*> common_subtrees;
};



#endif /* terracetree_hpp */
