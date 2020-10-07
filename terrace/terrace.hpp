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
    
    int taxa_num;
    int part_num;
    bool flag_trees_linked;
    
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
     *  TODO: Map "epsilon" branches, which were not mapped by the general tree linking function
     */
    void linkTreesEmptyImage();
    void linkTreeEmptyImage(int part, NodeVector &part_taxa, TerraceNode *node = nullptr, TerraceNode *dad = nullptr);

    /*
     *  Print Info about the branch between parent tree and induced partition trees
     */
    void printMapInfo();
    
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
    void getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, vector<TerraceNeighbor*> *nei1_vec, vector<TerraceNeighbor*> *nei2_vec);
    
    /*
     *  Insert a new taxon to the parent tree, update induced partition trees, update mapping
     */
    
    void insertNewTaxon(string node_name, TerraceNeighbor *nei_1, TerraceNeighbor *nei_2);
    
    /*
     *  Clean all link neighbours and taxa on parent tree and on induced partition trees
     */
    
    void cleanAllLinkNeighboursAndTaxa();
    
    /*
     *  Prepare top-low induced partition tree pairs: induced tree from the terrace and a common subtree with the initial tree (to be modified by inserting new taxa). Top level provided by the passed terrace, low level by the current terrace (which is initial terrace).
     */
    
    void create_Top_Low_Part_Tree_Pairs(vector<Terrace*> &part_tree_pairs, Terrace *terrace);
};

#endif /* terrace_hpp */
