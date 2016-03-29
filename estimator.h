//
//  estimator.h
//  iqtree
//
//  Created by Olga on 21/03/16.
//
//

#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <stdio.h>
#include "iqtree.h"
#include "alignment.h"

/**
    Improved estimation of pattern probabilities using shrinked estimator
 
	@author Olga Chernomor <olga.chernomor@univie.ac.at>
 */

void estimatorAnalysis(Params* params, Alignment* alignment, IQTree* tree);


/**
 * print pattern info (log likelihoods, frequency) to a file
 * @param filename output file name
 * @param tree phylogenetic tree
 * @param ptn_lh pattern log-likelihoods, will be computed if NULL
 * @param append TRUE to append to existing file, FALSE otherwise
 * @param linename name of the line, default "Site_Lh" if NULL
 */
void printPatternLhFreq(const char*filename, PhyloTree *tree, double *ptn_lh = NULL,
                 bool append = false, const char *linename = NULL);




#endif
