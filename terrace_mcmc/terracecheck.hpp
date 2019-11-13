//
//  terracecheck.hpp
//  iqtree
//
//  Created by Olga on 05.11.19.
//

#ifndef terracecheck_hpp
#define terracecheck_hpp

#include <stdio.h>
#include "utils/tools.h"
#include "tree/mtree.h"
#include "alignment/superalignment.h"

/*
    Check if a set of trees belongs to the terrace defined by a supermatrix and a user-input representative tree
 */

void runTerraceCheck(Params &params);


/*
 *  Get a corresponding induced partition tree
 */
MTree* getInducedPartitionTree(MTree* tree_complete, Alignment* aln, int part);

/*
 * Get a set of all induced partition trees
 */
void getALLInducedPartitionTrees(MTree* tree, Alignment* aln, vector<MTree*> &induced_part_trees);

/*
 * Get a tree topology as a string
 */
string getTreeTopologyString(MTree* tree);

/*
 *  Check if a tree lies on the terrace defined by a set of induced partition trees
 */

int terrace_check_two_trees(vector<MTree*> &induced_part_trees_rep, MTree* tree_query);

/*
 *  Check if all the taxa appear in the alignment
 */
bool checkTaxonSetsAlnTree(Alignment* aln, MTree* tree);

/*
 * Convert rooted tree to unrppted one. The same function is available for PhyloTree, but since the usage is here on MTree, the alignmnet, which is used in the original function, is not available
 */
void convertToUnrooted_MTree(MTree* tree,int nseq);

#endif /* terracecheck_hpp */
