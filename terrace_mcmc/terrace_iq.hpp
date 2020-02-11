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
    
    /**
     destructor
     */
    ~Terrace_IQ();
    
    int n_part;
    Alignment* aln;
    vector<MTree*> induced_part_trees;
    
    // TODO: while reading presence-absence matrix also fill out the vector of taxa_names, so the row ids correspond to correct taxa
    vector<IntVector> presence_absence_matrix;
    vector<string> taxa_names;
    
    
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
    
};
#endif /* terrace_iq_hpp */
