/*
 * node.h
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */

#pragma once

#include <string>
#include <set>
#include <map>
#include <bitset>
#include <vector>
#include "supermatrix.h"
// #include <boost/uuid/sha1.hpp>
#include <boost/uuid/detail/sha1.hpp>

const int MAX_TAXA = 1320;
const int MAX_LOCI = 1200;

class Supermatrix;

class MCNode {
public:
    MCNode();

    MCNode* addChild();

    void swapParents(MCNode* a, MCNode* b);

    bool isRoot();

    bool isLeaf();

    void debug(std::string indent = "");

    std::bitset<MAX_TAXA>* getLeaves();

    void getTaxonNames(std::set<std::string> &names);

    std::bitset<MAX_LOCI>* getLoci();

    void setName(std::string name);    

    MCNode* otherChild(MCNode* child);

    void replaceChild(MCNode* child, MCNode* replacement);

    void replaceChild2(MCNode* child, MCNode* replacement);

    void deleteChild(MCNode* child);

    void printNeighbours();

    void performNNI(MCNode* nni[4], bool variant, bool noreset = false);

    void root(MCNode* root, std::string taxon);

    void setTaxonIndices(std::map<std::string, int>* index);

    void setTaxonLoci(std::map<int, std::bitset<MAX_LOCI>*>* locus_by_taxon);

    void getOnTerraceNNIs(Supermatrix* supermatrix); 

    void getAllNNIs(); 

    void getNNI(int no, MCNode* nni[4]);

    std::string print(bool rooted = true, bool weights = false);

    std::string printSorted(bool rooted = true, bool weights = false);

    std::string getName();

    std::string name = "";

    std::string getFirstLeaf();

    int name_index;

    bool is_leaf = false;

    bool is_root = false;

    bool leaves_reset = true;

    void setRoot(bool root = true);

    std::string toString();

    bool hasCommonLeaves(std::bitset<MAX_TAXA>* a, bool complement = false, MCNode* n = NULL);

    MCNode* parent = NULL;

    MCNode* left_child = NULL;

    MCNode* right_child = NULL;   

    MCNode* clone();

    std::bitset<MAX_TAXA>* recalculateLeaves();

    std::bitset<MAX_LOCI>* recalculateLoci();

    void resetLeaves();

    std::string getTerraceId(Supermatrix* sm);

    void resetNNICache();

    void replaceLeafNames(std::map<int,int>);

    void getLeafNodes(std::vector<MCNode*> &nodes);

    std::string getSupertreeString(std::string supertreeString, std::map<int, std::map<int,int>>, std::map<int, std::map<int,int>>);

    std::string getSupertreeString(std::string supertreeString);

    std::string getSupertreeString(Supermatrix* sm);

    void getSubtree(std::set<std::string> &taxa);

    ~MCNode();

    MCNode* nni_cache[4*(MAX_TAXA-3)];
    int num_nnis = -1;

    void sort();  

protected:
    std::bitset<MAX_TAXA>* leaves = new std::bitset<MAX_TAXA>();

    //std::bitset<MAX_LOCI>* loci; //= new std::bitset<MAX_LOCI>();

    void _getOnTerraceNNIs(Supermatrix* supermatrix, int &result_count, MCNode* result[4*(MAX_TAXA-3)]);

    void _getAllNNIs(int &result_count, MCNode* result[4*(MAX_TAXA-3)]);
};
