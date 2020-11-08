//
//  presenceabsencematrix.hpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#ifndef presenceabsencematrix_hpp
#define presenceabsencematrix_hpp

#include <stdio.h>
#include "utils/tools.h"
#include "tree/mtree.h"
#include "terracetree.hpp"

class PresenceAbsenceMatrix {
    
public:
    
    vector<IntVector> pr_ab_matrix;
    vector<string> taxa_names;
    
    int taxa_num;
    int part_num;
    double missing_percent;
    
    /*
     *  initialization
     */
    void init();

    /*
     *  Reading a presence-absence matrix for supermatrix
     */
    void read_pr_ab_matrix(const char *infile);
    void read_pr_ab_matrix(istream &in);
    
    void print_pr_ab_matrix();
    
    /*
     * Find taxon id by taxon name
     */
    
    int findTaxonID(string taxon_name);
    
    /*
     *  a NodeVector part_taxa of (pointers to) leaf nodes from the tree is filled according to the presence-absence info for part
     */
    void getPartTaxa(int part, MTree *tree, MTree *part_tree, NodeVector &part_taxa);
    
    /*
     *  reorder vector of taxa_names and pr_ab_matrix according to the the order of taxa on a tree (based on node->ids)
     */
    void reorderAccordingToTree(NodeVector taxa_nodes);
    
    /*
     *   A variable to keep track if the reordering of taxa according to the tree was already performed
     */
    bool flag_reorderAccordingToTree;
    
    /*
     *  Get a submatrix of presence-absence matrix corresponding to taxa passed through taxa_names.
     */
    
    void getSubPrAbMatrix(vector<string> taxa_names_subset, PresenceAbsenceMatrix *submatrix, IntVector *parts = nullptr);
    void getSubPrAbMatrix(NodeVector taxon_nodes, PresenceAbsenceMatrix *submatrix, IntVector *parts = nullptr);
    
    /*
     *  Extend presence-absence matrix by the new taxon
     */
    
    void extend_by_new_taxa(string taxon_name, IntVector pr_ab_pattern);
    
    /*
     *  Remove a taxon from presence_absence matrix
     */
    
    void remove_taxon(string taxon_name);
    
    /*
     *  Function to get informaton about initial tree and a taxon order for stepwise insertion.
     *  taxa_names_sub - list of taxa to be on the intial tree (TODO: maybe get the tree right away?).
     *  list_taxa_to_insert - taxa in the order to be inserted
     */
    
    void getINFO_init_tree_taxon_order(vector<string> &taxa_names_sub, vector<string> &list_taxa_to_insert);
    
    /*
     *  Order taxa by their coverage
     */
    
    void orderTaxaByCoverage(vector<int> &taxon_ids, vector<IntVector> &coverage_info, IntVector &ordered_taxa);
    
    /*
     *  Compute amount of missing sequences (in %)
     */
    
    void percent_missing();
};

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree);

#endif /* presenceabsencematrix_hpp */
