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
#include "mtree.h"

/**
 *      main function to run clade analysis
 */
void runCladeAnalysis();

/**
 *      Analysis of Clades: Topology of the smallest clade subtree containing all species of interest
 */

class CladeAnalysis : public MTree {
    
public:
    /**
     *      Constructor
     */
     CladeAnalysis(MTree* tree);
    ~CladeAnalysis();
    
    /**
     *      Main function to start the analysis
     */

    void startCladeAnalysis(MTree* tree);
    
    /**
     *      Initialize variables
     */
    void initCladeAnalysis(MTree* tree);

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
    void checkClade(vector<int> *taxaSplit, bool *foundALL, bool *foundSOME);
    
    /**
     *      Set details about the smallest clade that contains all taxa from the list
     *      @param tree - input tree on which we perform the clade analysis
     *      @param taxaSplit - current smallest clade
     */
    void setMinClade(MTree *tree, vector<int> *taxaSplit);
    
    /**
     *      Print the results: 
     *      - smallest clade containing all taxa of interest (subtree)
     *      - its size
     *      - list of species on this clade
     *      - is it the smallest possible clade? (i.e. clade contains only species of interest)
     */
    void printResultsCA();
    
    /**
     *      Checks whether minCladeSize = taxaNameNUM, i.e. whether clade contains only species of interest
     */
    bool minPossible();

    
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
     *      The number of Species to analyse the clade for
     */
    int taxaNameNUM = NULL;
    
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
