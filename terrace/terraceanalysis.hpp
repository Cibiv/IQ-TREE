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

/*
 *  Main function for terrace analysis
 */
void runterraceanalysis(Params p);

/*
 *  Reading a presence-absence matrix for supermatrix
 */
IntVector read_pr_ab_matrix(const char *file);

#endif /* terraceanalysis_hpp */
