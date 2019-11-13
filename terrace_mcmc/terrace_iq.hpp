//
//  terrace_iq.hpp
//  iqtree
//
//  Created by Olga on 07.11.19.
//

#ifndef terrace_iq_hpp
#define terrace_iq_hpp

#include <stdio.h>
#include "tree/mtree.h"
#include "alignment/superalignment.h"

class Terrace_IQ {
    
public:
    /**
     constructor
     */
    Terrace_IQ();
    
    /**
     constructor
     */
    Terrace_IQ(Alignment* alignment, MTree* tree);
    
    /**
     constructor
     */
    Terrace_IQ(Alignment* alignment, vector<MTree*> trees);
    
    /**
     destructor
     */
    ~Terrace_IQ();
    
    int n_part;
    Alignment* aln;
    vector<MTree*> induced_part_trees;
    
    /*
     *  Modify terrace by setting a new set of induced partition trees
     */
    void setInducedPartitionTrees(vector<MTree*> trees);
    
    /*
     * Check if a tree is from the terrace
     */
    bool checkTree(MTree* tree);
    
};
#endif /* terrace_iq_hpp */
