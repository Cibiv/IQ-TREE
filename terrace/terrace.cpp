//
//  terrace.cpp
//  iqtree
//
//  Created by Olga on 14.09.20.
//

#include "terrace.hpp"
#include "terracenode.hpp"

Terrace::Terrace(){};
Terrace::~Terrace(){};

void Terrace::init(){
    
    taxa_num = 0;
    part_num = 0;
    flag_trees_linked = false;
}

Terrace::Terrace(const char *infile_tree, bool is_rooted,const char *infile_matrix){
    
    init();
    
    readTree(infile_tree,is_rooted);
    
    matrix = new PresenceAbsenceMatrix();
    matrix->read_pr_ab_matrix(infile_matrix);
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();
    
    get_part_trees();
    printInfo();
};

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix m){
    
    init();

    copyTree(&tree);
    matrix = &m;
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();
    
    get_part_trees();
    printInfo();
}

void Terrace::get_part_trees(){
    
    assert(matrix!=nullptr && "ERROR: Presence-absence matrix is absent! I cannot get induced partition trees..");
    cout<<"Preparing induced partition trees..."<<endl;
    
    int id;
    string taxa_set = "";
    string taxon_name;
    NodeVector taxa_nodes;
    
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(taxa_num);
    
    getTaxa(taxa_nodes);
    for(int part=0; part<part_num; part++){
        TerraceTree* induced_tree = new TerraceTree();
        for(it2=taxa_nodes.begin();it2!=taxa_nodes.end();it2++){
            taxon_name=(*it2)->name;
            id=matrix->findTaxonID(taxon_name);
            assert(id != -1 && "Not all of the taxa appear in the pr_ab_matrix!");
            check_int[id] = matrix->pr_ab_matrix[id][part];
            //cout<<"Taxon["<<id<<"] in partition "<<part<<" is "<<check_int[id]<<endl;
        }
        taxa_set.clear();
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        induced_tree->copyTree(this,taxa_set);
        induced_trees.push_back(induced_tree);
        //induced_trees[part]->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    }
    
}

void Terrace::printInfo(){

    cout<<endl<<"Printing terrace information:"<<endl<<endl;
    cout<<"Terrace representative tree:"<<endl;
    this->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    cout<<endl;
    
    matrix->print_pr_ab_matrix();
    
    if(induced_trees.size()>0){
        int i=0;
        cout<<"Induced partition trees:"<<endl;
        for(vector<TerraceTree*>::iterator it = induced_trees.begin(); it < induced_trees.end(); it++){
            i++;
            cout<<"Part["<<i<<"]: ";
            (*it)->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
            cout<<endl;
        }
    }
}

void Terrace::linkTrees(){
    NodeVector part_taxa;

    if(!matrix->flag_reorderAccordingToTree){
        NodeVector tree_taxa;
        this->getTaxa(tree_taxa);
        matrix->reorderAccordingToTree(tree_taxa);
    }
    
    // Map only branches from a parent tree to induced partition trees, equivalent to PTA
    for(int part=0; part<part_num; part++){
        part_taxa.clear();
        matrix->getPartTaxa(part, this, induced_trees[part], part_taxa);
        
        // TODO:
        // - you need to make sure that you only map to induced partition trees if they have more than two leaves in common with parent tree
        // - if the number of leaves is equal or less than two, than the tree does not provide any constraint for the extention of the tree by new taxa
        
        linkTree(part,part_taxa);
    }
    flag_trees_linked = true;
    
}

void Terrace::linkTree(int part, NodeVector &part_taxa, TerraceNode *node, TerraceNode *dad){
    
    if (!node) {
        if (!root->isLeaf())
            node = (TerraceNode*) root;
        else
            node = (TerraceNode*) root->neighbors[0]->node;
        ASSERT(node);
        if (node->isLeaf()) // two-taxa parent tree
            dad = (TerraceNode*)node->neighbors[0]->node;
    }
    TerraceNeighbor *nei = NULL;
    TerraceNeighbor *dad_nei = NULL;
    if (dad) {
        nei = (TerraceNeighbor*)node->findNeighbor(dad);
        dad_nei = (TerraceNeighbor*)dad->findNeighbor(node);
        if (nei->link_neighbors.empty()) nei->link_neighbors.resize(part_num);
        if (dad_nei->link_neighbors.empty()) dad_nei->link_neighbors.resize(part_num);
        nei->link_neighbors[part] = NULL;
        dad_nei->link_neighbors[part] = NULL;
    }
    if (node->isLeaf()) {
        ASSERT(dad);
        TerraceNode *node_part = (TerraceNode*)part_taxa[node->id]; //part_taxa[] is a node on PARTITION tree
        if (node_part) {
            TerraceNode *dad_part = (TerraceNode*)node_part->neighbors[0]->node;
            TerraceNeighbor *dad_part_nei = (TerraceNeighbor*)dad_part->findNeighbor(node_part);
            ASSERT(node_part->isLeaf());
            nei->link_neighbors[part] = (TerraceNeighbor*) node_part->neighbors[0];
            dad_nei->link_neighbors[part] = dad_part_nei;
            
            // TODOOOOOOOOOOO
            // Maybe a separate function for the back synchronisation would be better, but check carefully.
            // Actually, for generation you need kind of a different map.
            // You only need a map from induced partition trees to parent tree and, in principle, you should not care about the other map...
            
            // Use this one function to do all the maps you need
            // Back map from induced partition tree onto parent tree
            // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
            //((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors.push_back(nei);
            //dad_part_nei->link_neighbors.push_back(dad_nei);
        } else {
            
            // the empty branches are needed for a map from parent to partition
            dad->empty_branches.push_back(node->neighbors[0]->id);
            dad->empty_br_dad_nei.push_back(dad->findNeighbor(node));
            dad->empty_br_node_nei.push_back(node->findNeighbor(dad));
            
            // the empty taxa are needed for a map from induced partition tree to "common" induced partition tree
            dad->empty_taxa.push_back(node);
        }
        return;
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        linkTree(part, part_taxa, (TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
    if (!dad) return;
    
    linkBranch(part, nei, dad_nei);
}

void Terrace::linkBranch(int part, TerraceNeighbor *nei, TerraceNeighbor *dad_nei) {
    TerraceNode *node = (TerraceNode*)dad_nei->node;
    TerraceNode *dad = (TerraceNode*)nei->node;
    nei->link_neighbors[part] = NULL;
    dad_nei->link_neighbors[part] = NULL;
    vector<TerraceNeighbor*> part_vec;
    vector<TerraceNeighbor*> child_part_vec;
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        if (((TerraceNeighbor*)*it)->link_neighbors[part]) {
            part_vec.push_back((TerraceNeighbor*)(((TerraceNeighbor*)*it)->link_neighbors[part]));
            child_part_vec.push_back((TerraceNeighbor*)(((TerraceNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]));
            ASSERT(child_part_vec.back()->node == child_part_vec.front()->node || child_part_vec.back()->id == child_part_vec.front()->id);
        }
    }
    
    // node or dad empty branch/taxa?
    int node_empty_branch_size = node->empty_branches.size();
    int node_empty_taxa_size = node->empty_taxa.size();
    
    if (part_vec.empty()){
        if(node_empty_branch_size>0){
            for(int i=0; i<node_empty_branch_size; i++){
                dad->empty_branches.push_back(node->empty_branches[i]);
                dad->empty_br_dad_nei.push_back(node->empty_br_dad_nei[i]);
                dad->empty_br_node_nei.push_back(node->empty_br_node_nei[i]);
            }
            //node->empty_branches.clear();
        }
        dad->empty_branches.push_back(nei->id);
        
        if(node_empty_taxa_size>0){
            for(int i=0; i<node_empty_taxa_size; i++){
                dad->empty_taxa.push_back(node->empty_taxa[i]);
            }
            //node->empty_taxa.clear();
        }
        return;
    }
    
    int i=0;
    if (part_vec.size() == 1) {
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = part_vec[0];
        
        if(node_empty_branch_size>0){
            for(i=0; i<node_empty_branch_size; i++){
                ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part]=child_part_vec[0];
                ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = part_vec[0];
                
                // TODO: actually, you need a map back: for child_part_vec[0] and part_vec[0] in linked_neighbors save all branche sthat map to them
            }
        }
        
        if(node_empty_taxa_size>0){
            for(i=0; i<node_empty_taxa_size; i++){
                child_part_vec[0]->taxa_to_insert.push_back(node->empty_taxa[i]);
                part_vec[0]->taxa_to_insert.push_back(node->empty_taxa[i]);
            }
            node->empty_taxa.clear();
        }
        
        if(dad->empty_taxa.size()>0){
            for(i=0; i<dad->empty_taxa.size(); i++){
                child_part_vec[0]->taxa_to_insert.push_back(dad->empty_taxa[i]);
                part_vec[0]->taxa_to_insert.push_back(dad->empty_taxa[i]);
            }
            dad->empty_taxa.clear();
        }
        
        return;
    }
    
    // what's this below? =´(
    if (part_vec[0] == child_part_vec[1]) {
        // ping-pong, out of sub-tree
        ASSERT(part_vec[1] == child_part_vec[0]);
        return;
    }
    
    TerraceNode *node_part = (TerraceNode*) child_part_vec[0]->node;
    TerraceNode *dad_part = NULL;
    FOR_NEIGHBOR(node_part, NULL, it) {
        bool appear = false;
        for (vector<TerraceNeighbor*>::iterator it2 = part_vec.begin(); it2 != part_vec.end(); it2++){
            if ((*it2) == (*it)) {
                appear = true; break;
            }
        }
        if (!appear) {
            ASSERT(!dad_part);
            dad_part = (TerraceNode*)(*it)->node;
        }
    }
    nei->link_neighbors[part] = (TerraceNeighbor*)node_part->findNeighbor(dad_part);
    dad_nei->link_neighbors[part] = (TerraceNeighbor*)dad_part->findNeighbor(node_part);
}

void Terrace::printMapInfo(){
    
    cout<<"Mapping info"<<endl<<endl;
    NodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    int part = 0;
    drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    for (vector<TerraceTree*>::iterator it = induced_trees.begin(); it != induced_trees.end(); it++, part++) {
        cout << "Subtree for partition " << part << endl;
        (*it)->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE );
        for (int i = 0; i < nodes1.size(); i++) {
            Neighbor *nei1 = ((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part];
            Neighbor *nei2 = ((TerraceNeighbor*)nodes2[i]->findNeighbor(nodes1[i]))->link_neighbors[part];
            cout << nodes1[i]->findNeighbor(nodes2[i])->id << ":";
            if (nodes1[i]->isLeaf()) cout << nodes1[i]->name; else cout << nodes1[i]->id;
            cout << ",";
            if (nodes2[i]->isLeaf()) cout << nodes2[i]->name; else cout << nodes2[i]->id;
            cout <<"("<<nodes1[i]->findNeighbor(nodes2[i])->length<<")"<< " -> ";
            if (nei2) {
                cout << nei2->id << ":";
                if (nei2->node->isLeaf())
                    cout << nei2->node->name;
                else cout << nei2->node->id;
            }
            else cout << -1;
            cout << ",";
            if (nei1){
                if (nei1->node->isLeaf())
                    cout << nei1->node->name;
                else cout << nei1->node->id;
                cout <<"("<<nei1->length<<")";
            }
            else cout << -1;
            cout << endl;
        }
    }
}


void Terrace::linkTreesEmptyImage(){
    
    assert(flag_trees_linked && "The parent and the induced partition trees are not linked yet! Exit.");
    NodeVector part_taxa;
    
    for(int part=0; part<part_num; part++){
        part_taxa.clear();
        matrix->getPartTaxa(part, this, induced_trees[part], part_taxa);
        linkTreeEmptyImage(part,part_taxa);
    }
}

void Terrace::linkTreeEmptyImage(int part, NodeVector &part_taxa, TerraceNode *node, TerraceNode *dad){
    
};



