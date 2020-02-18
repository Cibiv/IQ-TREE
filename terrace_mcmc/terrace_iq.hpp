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
     *  TODO: add another constructor with presence_abscence_matrix and modify accordingly parameters and functions that use functionality of aln
     */
    
    Terrace_IQ(const char* file_presence_absence, MTree* tree);
    Terrace_IQ(vector<IntVector> matrix, vector<string> names, MTree* tree);
    
    /**
     destructor
     */
    ~Terrace_IQ();
    
    int n_part;
    int taxa_num;
    Alignment* aln;
    MTree* representative_tree;
    vector<MTree*> induced_part_trees;
    
    // TODO: while reading presence-absence matrix also fill out the vector of taxa_names, so the row ids correspond to correct taxa
    vector<IntVector> pr_ab_matrix;
    vector<string> taxa_names;
    
    /*
     *  If terrace is constructed from an alignment instead of presence-absence matrix, initialize taxa_names and pr_ab_matrix
     */
    void initTaxaNamesPresenceAbsenceM();
    
    
    /*
     *  Modify terrace by setting a new set of induced partition trees
     */
    void setInducedPartitionTrees(vector<MTree*> trees);
    
    /*
     * Check if a tree is from the terrace
     */
    bool checkTree(MTree* tree);
    
    /*
     * Print out induced partition trees
     */
    void printInducedTrees(ostream &out);
    
    /*
     * Read presence-absence matrix from a file
     */
    void readPresenceAbsenceMatrix(const char *infile);
    void readPresenceAbsenceMatrix(istream &in);
    
    /*
     * Find a taxon ID in pr_ab_matrix
     */
    int getTaxonID_in_pr_ab_m(string taxon_name);
    
    /*
     *  Get induced partition trees using pr_ab_matrix
     */
    void getALLInducedPartitionTreesM();
    
};
#endif /* terrace_iq_hpp */
