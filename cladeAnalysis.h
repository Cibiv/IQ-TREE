//
//  cladeAnalysis.h
//  iqtree
//
//  Created by Olga on 31/08/16.
//
//

#ifndef __iqtree__cladeAnalysis__
#define __iqtree__cladeAnalysis__

#include <stdio.h>
#include "iqtree.h"

/**
 *      Analysis of Clades: Topology of the smallest clade subtree containing all species of interest
 */

class CladeAnalysis : public MTree {
    
public:
    /**
     *      Constructor
     */
     CladeAnalysis(IQTree* tree);
    ~CladeAnalysis();
    
    /**
     *      Main function to start the analysis
     */

    void startCladeAnalysis(IQTree* tree);
    
    /**
     *      Initialize variables
     */
    void initCladeAnalysis(IQTree* tree);

    /**
     *      Reading input list of taxa to analyse the clade for
     */
    void readInputTaxa(const char* infile);
    
    /**
     *      Reading input list of taxa to analyse the clade for
     */
    void readInputTaxa(istream &in);
    
    /**
     *      Function to check whether all taxa of interest are present in the current set A (split A|B) of species
     *      @param taxaSplit - set of taxa from one side of the current split A|B, e.g. A
     *      @param foundALL - bool variable to track whether all taxa of interest were found in subset A
     *      @param foundSOME - bool variable to track whether at least some of the taxa of interest were found in A
     */
    void checkClade(vector<string> *taxaSplit, bool *foundALL, bool *foundSOME);
    
    /**
     *      Set details about the smallest clade that contains all taxa from the list
     *      @param tree - input tree on which we perform the clade analysis
     *      @param taxaSplit - current smallest clade
     */
    void setMinClade(IQTree *tree, vector<string> *taxaSplit);

    
private:
    
    /**
     *      Names of Species to analyse the clade for
     */
    vector<string> taxaName;
    
    /**
     *      IDs of Species to analyse the clade for
     */
    vector<int> taxaNameID;
    
    /**
     *      The size of the smallest clade containing all species of interest
     */
    
    int minCladeSize;
    
    /**
     *      List of Species of the smallest clade containg all species of interest
     */
    
    vector<string> minCladeSpecies;
    
};

#endif /* defined(__iqtree__cladeAnalysis__) */
