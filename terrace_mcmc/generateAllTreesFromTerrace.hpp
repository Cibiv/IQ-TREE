//
//  generateAllTreesFromTerrace.hpp
//  iqtree
//
//  Created by Olga on 10.02.20.
//

#ifndef generateAllTreeFromTerrace_hpp
#define generateAllTreeFromTerrace_hpp

#include <stdio.h>
#include "terrace_iq.hpp"
#include "terracecheck.hpp"

#endif /* generateAllTreeFromTerrace_hpp */

/*
 *  Here is the algorithm to generate all trees from a terrace for a general case (a reference taxon is not required).
 */

void generateTerraceTrees(Terrace_IQ* terrace, const char *userOutFile);
