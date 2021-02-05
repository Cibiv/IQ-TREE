//
//  terrace.hpp
//  iqtree
//
//  Created by Olga on 14.09.20.
//

#ifndef terrace_hpp
#define terrace_hpp

#include <stdio.h>
#include "terracetree.hpp"
#include "terracenode.hpp"
#include "presenceabsencematrix.hpp"

class Terrace: public TerraceTree
{
public:
    /*
     * constructor
     */
    Terrace();
    
    /*
     * constructor
     */
    Terrace(TerraceTree tree, PresenceAbsenceMatrix *m);
    
    /*
     * constructor
     */
    Terrace(TerraceTree tree, PresenceAbsenceMatrix *m, vector<TerraceTree*> input_induced_trees);
    
    /*
     *  constructor
     */
    Terrace(const char *infile_tree, bool is_rooted,const char *infile_matrix);
    
    /*
     * destructor
     */
    ~Terrace();
    
    /*
     *  Initialisator
     */
    void init();
    
    vector<TerraceTree*> induced_trees;
    PresenceAbsenceMatrix *matrix;
    StrVector terrace_trees;
    
    int taxa_num;
    int part_num;
    
    /*
     * Original taxon names
     */
    StrVector taxa_names_orgn;
    
    /*
     *  Number of trees on terrace
     */
    int terrace_trees_num;
    
    /*
     *  Number of intermediate trees visited
     */
    int intermediated_trees_num;
    
    /*
     *  Number of dead ends encountered
     */
    int dead_ends_num;
    
    // file to output all generated terrace trees
    string out_file;
    bool terrace_out;
    
    // Stopping rules
    int terrace_max_trees;
    int intermediate_max_trees;
    int seconds_max;
    
    /*
     *  Print terrace info: a representative tree, induced trees and presence-absence matrix
     */

    void printInfo();
    
    /*
     *  get induced partition trees
     */
    void get_part_trees();
    
    /*
     *  set induced partition trees, if they already exist
     */
    void set_part_trees(vector<TerraceTree*> input_induced_trees);
    
    /*
     * link parent tree and induced partition trees
     */
    
    void linkTrees(bool back_branch_map, bool back_taxon_map);
    
    /*
     * link parent tree and a single induced partition tree
     */
    void linkTree(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node = nullptr, TerraceNode *dad = nullptr);
    
    /*
     *  link one branch of parent tree on partition tree part (for internal branches only, oder?)
     */
    void linkBranch(int part, TerraceNeighbor *nei, TerraceNeighbor *dad_nei, bool back_branch_map, bool back_taxon_map);
    
    /*
     *  Local map update after insertion/deletion of the taxon
     */
    
    void update_map(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad = nullptr);
    
    /*
     *  Print Info about the branch between parent tree and induced partition trees
     */
    void printMapInfo(int partition = -1);
    
    /*
     *  Print Info about inverse branch and taxon images from induced partition trees to master and upper level induced trees, respectivelly
     */
    void printBackMapInfo();
    
    /*
     *  For all nodes of the tree clear information about empty branches and empty nodes, which might have remained after the previous partition
     */
    
    void clearEmptyBranchAndTaxaINFO(TerraceNode *node = nullptr, TerraceNode *dad = nullptr);
    
    /*
     *  For a given taxon name get allowed branches. aux_terrace contains pointers to top level induced partition trees (which are stored as terraces)
     */
    void getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, NodeVector *node1_vec, NodeVector *node2_vec);
    
    /*
     *  Insert a new taxon to the parent tree, update induced partition trees, update mapping (locally)
     */
    void extendNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> aux_terrace,vector<int> pr_ab_info);
    
    /*
     *  Insert a new taxon to the parent tree, update induced partition trees, update mapping. The maps are not updated, but cleaned and then everything is relinked again from scratch
     */
    void extendNewTaxon_naive(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> part_tree_pairs,vector<int> pr_ab_info);
    
    /*
     *  Clean all link neighbours and taxa on parent tree and on induced partition trees
     */
    
    void cleanAllLinkNeighboursAndTaxa(bool clean_induced_part_maps = false);
    
    /*
     *  Prepare top-low induced partition tree pairs: induced tree from the terrace and a common subtree with the initial tree (to be modified by inserting new taxa). Top level provided by the passed terrace, low level by the current terrace (which is initial terrace).
     */
    
    void create_Top_Low_Part_Tree_Pairs(vector<Terrace*> &part_tree_pairs, Terrace *terrace);
    
    /*
     *  The main function to generate trees by recursive taxon insertion
     */
    
    void generateTerraceTrees(Terrace *terrace, vector<Terrace*> part_tree_pairs, vector<string> *list_taxa_to_insert, int taxon_to_insert, bool *progress_status);
    
    /*
     *  Remove one taxon from the terrace tree, from induced partition trees, update the mapping
     *  Here the update is local update: only involved branches and maps of trees are beeing updated.
     */
    
    void remove_one_taxon(string taxon_name, vector<Terrace*> part_tree_pairs);
    
    /*
     *  Remove one taxon from the terrace tree, from induced partition trees, update the mapping.
     *  Here the update corresponds to:
     *  - remove taxon,
     *  - clear all maps from agile (parent) tree and from low-top induced partition trees,
     *  - relink all branches of all trees from scratch
     */
    void remove_one_taxon_naive(string taxon_name, vector<Terrace*> part_tree_pairs);
    
    /*
     *  Re-link
     */
    void relinkALL(vector<Terrace*> part_tree_pairs);
    
    /*
     * Print initial tree, print TOP and LOW induced trees
     */
    void print_ALL_DATA(vector<Terrace*> part_tree_pairs);
    
    /*
     *  Rename taxa on a tree and in presence-absence matrix
     */
    void renameTaxa();
    
    /*
     *  Check if the query tree has the same set of induced partition trees
     */
    bool check_two_trees(MTree* query_tree);
    
    /*
     *  Write all generated trees to file
     */
    void write_terrace_trees_to_file();
    
    /*
     *  Write summary of generating trees from a terrace
     */
    void write_summary_generation();
    
    /*
     *  Write warning depending on the activated stopping rule
     */
    void write_warning_stop(int type);
};

#endif /* terrace_hpp */
