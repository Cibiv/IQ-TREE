//
//  terracetree.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "terracetree.hpp"
#include "terracenode.hpp"

TerraceTree::TerraceTree(): MTree(){};
TerraceTree::~TerraceTree(){};


void TerraceTree::readTree(const char *infile, bool &is_rooted){
    MTree::readTree(infile, is_rooted);
};

void TerraceTree::readTree(istream &in, bool &is_rooted) {
    MTree::readTree(in, is_rooted);
}

Node* TerraceTree::newNode(int node_id, const char* node_name) {
    return (Node*) (new TerraceNode(node_id, node_name));
}

Node* TerraceTree::newNode(int node_id, int node_name) {
    return (Node*) (new TerraceNode(node_id, node_name));
}


void TerraceTree::copyTree_byTaxonNames(MTree *tree, vector<string> taxon_names){
    
    //cout<<"Copying a tree using a vector of taxon names..."<<endl;
    
    int i,j, sum = 0;
    int taxa_num = taxon_names.size();
    bool check;
    
    string taxa_set = "";
    NodeVector taxa_nodes;
    
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(tree->leafNum);
    std::fill(check_int.begin(), check_int.end(), 0);
    
    tree->getTaxa(taxa_nodes);
    
    for(j=0; j<taxa_nodes.size(); j++){
        check = false;
        for(i=0; i<taxa_num; i++){
            if(taxa_nodes[j]->name == taxon_names[i]){
                check_int[j]=1;
                sum+=1;
                check = true;
                break;
            }
        }
        //cout<<"Taxon["<<j<<"] = "<<check_int[j]<<"| main_tree:"<<taxa_nodes[j]->name<<"-> subtree: "<<((check) ? taxon_names[i] : "")<<endl;
    }
    assert(sum == taxa_num && "Not all of the taxa appear in the complete tree!");
    taxa_set.clear();
    taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
    copyTree(tree,taxa_set);

    //printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
}

void TerraceTree::cleanAllLinkINFO(TerraceNode *node, TerraceNode *dad){
    
    if(!node){
        if(root->isLeaf()){
            node = (TerraceNode*) root->neighbors[0]->node;
        }else{
            node = (TerraceNode*) root;
        }
        
        ASSERT(node);
    }
    
    int part;
    if(dad){
        TerraceNeighbor *nei = (TerraceNeighbor*)node->findNeighbor(dad);
        TerraceNeighbor *dad_nei = (TerraceNeighbor*)dad->findNeighbor(node);
        if(nei->link_neighbors.size()>0){
            for(part=0; part<nei->link_neighbors.size(); part++){
                if(((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.size()>0){
                    ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.clear();
                    ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors.clear();
                }
                if(((TerraceNeighbor*)nei->link_neighbors[part])->taxa_to_insert.size()>0){
                    ((TerraceNeighbor*)nei->link_neighbors[part])->taxa_to_insert.clear();
                    ((TerraceNeighbor*)dad_nei->link_neighbors[part])->taxa_to_insert.clear();
                }
            }
            nei->link_neighbors.clear();
            dad_nei->link_neighbors.clear();
            
        }
        
        node->empty_taxa.clear();
        node->empty_branches.clear();
        node->empty_br_dad_nei.clear();
        node->empty_br_node_nei.clear();
        
        dad->empty_taxa.clear();
        dad->empty_branches.clear();
        dad->empty_br_dad_nei.clear();
        dad->empty_br_node_nei.clear();
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it){
        cleanAllLinkINFO((TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
}

void TerraceTree::insertNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch){
    
    // IMPORTANT!!!!!! WARNING: CAREFULL about calling a presence_absence matrix entries by the leaf id!!!!!! CHECK the corresponding code, because your maps will be fucked up!
    // TODO: safe option is a function to re-number taxa, nodes
    
    TerraceNode *node_1, *node_2;
    
    node_1 = (TerraceNode*) newNode(nodeNum, node_name.c_str());
    
    leafNum += 1;
    nodeNum += 1;
    
    node_2 = (TerraceNode*)(newNode(nodeNum));
    nodeNum += 1;
    
    int br_id = branchNum;
    branchNum += 1;
    
    // WARNING: you need unique branch ids for intersection in allowed branches!!!!!
    // CHECK: if there are any issues with current handling of branch ids
    node_1->addNeighbor(node_2, 0.0, br_id);
    node_2->addNeighbor(node_1, 0.0, br_id);
    
    br_id = node_1_branch->findNeighbor(node_2_branch)->id;
    node_1_branch->updateNeighbor(node_2_branch, node_2, 0.0);
    node_1_branch->findNeighbor(node_2)->id = br_id;
    node_2->addNeighbor(node_1_branch, 0.0, br_id);
    
    br_id = branchNum;
    node_2_branch->updateNeighbor(node_1_branch, node_2, 0.0);
    node_2_branch->findNeighbor(node_2)->id = br_id;
    
    node_2->addNeighbor(node_2_branch, 0.0, br_id);
    
    drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    
}
