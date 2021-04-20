//
//  terraceanalysis.hpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#ifndef terraceanalysis_hpp
#define terraceanalysis_hpp

#include <stdio.h>
#include "utils/tools.h"
#include "terrace/terrace.hpp"

/*
 *  Main function for terrace analysis
 */
void runterraceanalysis(Params &params);

/*
 * Run terrace check:
 * - check if query trees lie on the same terrace with a representative tree. Naive pairwise comparison of induced partition trees.
 */
void run_terrace_check(Terrace *terrace, Params &params);


/**
        Check if trees in the list have degree 2 nodes, if yes, collapse them. This is used to deal with trees output by terraphast.
 */
void run_collapse_trees(Params &params);

#endif /* terraceanalysis_hpp */
