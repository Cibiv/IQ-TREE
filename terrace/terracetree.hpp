//
//  terracetree.hpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#ifndef terracetree_hpp
#define terracetree_hpp

#include <stdio.h>
#include "mtree.h"

struct presence_absence_matrix{
    vector<string> taxa_names;
    vector<IntVector> pr_ab_matrix;
};


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

    
};



#endif /* terracetree_hpp */
