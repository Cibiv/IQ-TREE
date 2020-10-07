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
    
    NodeVector taxa_nodes;
    getTaxa(taxa_nodes);
    matrix->reorderAccordingToTree(taxa_nodes);
    
    get_part_trees();
    printInfo();
};

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix *m){
    
    init();

    copyTree(&tree);
    matrix = m;
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();
    
    NodeVector taxa_nodes;
    getTaxa(taxa_nodes);
    //printTree(cout);
    
    matrix->reorderAccordingToTree(taxa_nodes);
    //matrix->print_pr_ab_matrix();
    
    get_part_trees();
    printInfo();
}

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix *m, vector<TerraceTree*> input_induced_trees){
 
    init();
    
    copyTree(&tree);
    matrix = m;
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();
    
    NodeVector taxa_nodes;
    getTaxa(taxa_nodes);
    //printTree(cout);
    
    matrix->reorderAccordingToTree(taxa_nodes);
    //matrix->print_pr_ab_matrix();
    
    set_part_trees(input_induced_trees);
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

void Terrace::set_part_trees(vector<TerraceTree*> input_induced_trees){
    
    for(int i=0; i<input_induced_trees.size(); i++){
        induced_trees.push_back(input_induced_trees[i]);
    }
    
    // TODO: you need a check, that your input_induced_trees indeed correspond to presence-absence pattern of pr_ab_matrix
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

void Terrace::linkTrees(bool back_branch_map, bool back_taxon_map){
    NodeVector part_taxa;

    if(!matrix->flag_reorderAccordingToTree){
        NodeVector tree_taxa;
        this->getTaxa(tree_taxa);
        matrix->reorderAccordingToTree(tree_taxa);
    }
    
    for(int part=0; part<part_num; part++){
        part_taxa.clear();
        matrix->getPartTaxa(part, this, induced_trees[part], part_taxa);
        
        // TODO!!!!!!!!!!:
        // - you need to make sure that you only map to induced partition trees if they have more than two leaves in common with parent tree
        // - if the number of leaves is equal or less than two, than the tree does not provide any constraint for the extention of the tree by new taxa
        //   -> instead of calling linkTree(part, part_taxa) have a separate function for "no or less than two taxa" trees
        
        //cout<<"YUPEE-YUPEE-YEAH! PARTITION "<<part<<endl<<endl;
        //clearEmptyBranchAndTaxaINFO();
        linkTree(part, part_taxa, back_branch_map, back_taxon_map);
    }
    flag_trees_linked = true;
    
}

void Terrace::linkTree(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad){
    
    // SEHR WICHTIG! WARNING: do not mix mapping from the parent tree and upper level induced partition trees, because empty branches and empty taxa will be messed up on partition trees.
    
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
            
            if(back_branch_map){
                // Back map from induced partition tree onto parent tree
                // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors.push_back(nei);
                dad_part_nei->link_neighbors.push_back(dad_nei);
            }
            
        } else {
            
            // the empty branches are needed for a map from parent to partition
            dad->empty_branches.push_back(nei->id);
            dad->empty_br_dad_nei.push_back(dad_nei);
            dad->empty_br_node_nei.push_back(nei);
            
            if(back_taxon_map){
                // the empty taxa are needed for a map from induced partition tree to "common" induced partition tree
                // when you know the link_neighbor, you have to map all empty taxa to that neighbor
                dad->empty_taxa.push_back(node);
                //cout<<"INFO CHECK: the length of the empty leaves: dad("<<dad->id<<") = "<<dad->empty_taxa.size()<<endl;
            }
        }
        return;
    }
    
    //cout<<"MASTER NODE:"<<node->id<<endl;
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<endl;
        linkTree(part, part_taxa, back_branch_map, back_taxon_map, (TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
    if (!dad) {
        // Check, if the final dad has empty_branches and empty_taxa and map them to branch, which is available.
        // Note, that if there are some empty branches/taxa, there will be exactly one branch available for mapping.
        if(node->empty_taxa.size()>0 or node->empty_branches.size()>0){
            //cout<<"I got here! Partition: "<<part<<"| Node->id:"<<node->id<<endl;
            FOR_NEIGHBOR_DECLARE(node, NULL, it) {
                if(((TerraceNeighbor*)(*it))->link_neighbors[part]){
                    TerraceNeighbor* node_nei_part = (TerraceNeighbor*)((TerraceNeighbor*)(*it))->link_neighbors[part];
                    TerraceNeighbor* dad_nei_part = (TerraceNeighbor*)((TerraceNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part];
                    int i;
                    for(i=0; i<node->empty_br_dad_nei.size(); i++){
                        ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part] = node_nei_part;
                        ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = dad_nei_part;
            
                        if(back_branch_map){
                            node_nei_part->link_neighbors.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors.push_back(node->empty_br_dad_nei[i]);
                        }
                    }
                    
                    if(back_taxon_map){
                        //cout<<endl<<"INFO CHECK (in IF!dad): node="<<node->id<<" | empty taxa:"<<endl;
                        for(i=0; i<node->empty_taxa.size(); i++){
                            //cout<<" "<<node->empty_taxa[i]->id<<",";
                            node_nei_part->taxa_to_insert.push_back(node->empty_taxa[i]);
                            dad_nei_part->taxa_to_insert.push_back(node->empty_taxa[i]);
                        }
                        //cout<<endl;
                        node->empty_taxa.clear();
                    }
                        
                    node->empty_branches.clear();
                    node->empty_br_node_nei.clear();
                    node->empty_br_dad_nei.clear();
                    
                    return;
                }
            }
            
            // WARNING: I think, it should never rich this point, because at least one neighbour shuold have a corresponding branch on partition tree.
            cout<<"ERROR: hey, the assumption was, that the code should never come here..."<<endl;
        }
        return;
    }
    linkBranch(part, nei, dad_nei, back_branch_map, back_taxon_map);
}

void Terrace::linkBranch(int part, TerraceNeighbor *nei, TerraceNeighbor *dad_nei, bool back_branch_map, bool back_taxon_map) {
    //cout<<endl;
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
    
    // QUESTION: node or dad empty branch/taxa?
    // QUESTION: if you start mapping the dad, you might encounter situations, when a dad does not know yet, that there are some empty branches.
    // So at the moment you need only empty branches of node and the dad will be mapped at the final state (i.e. final node - the start of the tree traversal)
    
    int node_empty_branch_size = node->empty_branches.size();
    int node_empty_taxa_size = node->empty_taxa.size();
    
    if (part_vec.empty()){
        int i=0;
        //cout<<"CASE 0: CHILDREN DO NOT HAVE IMAGE"<<endl<<endl;
        if(node_empty_branch_size>0){
            for(int i=0; i<node_empty_branch_size; i++){
                dad->empty_branches.push_back(node->empty_branches[i]);
                dad->empty_br_dad_nei.push_back(node->empty_br_dad_nei[i]);
                dad->empty_br_node_nei.push_back(node->empty_br_node_nei[i]);
            }
            node->empty_branches.clear();
            node->empty_br_node_nei.clear();
            node->empty_br_dad_nei.clear();
        }
        
        // since below node the subtrees are empty, add also current branch to an empty list of the dad
        dad->empty_branches.push_back(nei->id);
        dad->empty_br_dad_nei.push_back(dad_nei);
        dad->empty_br_node_nei.push_back(nei);
        
        if(back_taxon_map){
            if(node_empty_taxa_size>0){
                /*cout<<"INFO CHECK (in two children nodes are empty): passing empty taxa to dad="<<dad->id<<" :"<<endl;
                if(dad->empty_taxa.size()>0){
                    cout<<"currently list of empty taxa for dad:";
                    for(i=0; i<dad->empty_taxa.size(); i++){
                        cout<<" "<<dad->empty_taxa[i]->id;
                    }
                    cout<<endl;
                }
                 */
                
                for(i=0; i<node_empty_taxa_size; i++){
                    //cout<<" "<<node->empty_taxa[i]->id;
                    dad->empty_taxa.push_back(node->empty_taxa[i]);
                }
                //cout<<endl;
                node->empty_taxa.clear();
            }
        }
        return;
    }
    
    int i=0;
    if (part_vec.size() == 1) {
        
        //cout<<"CASE 1: ONE CHILD HAS IMAGE"<<endl<<endl;
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = part_vec[0];
        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            part_vec[0]->link_neighbors.push_back(dad_nei);
        }
        
        // Get maps for empty branches
        if(node_empty_branch_size>0){
            //cout<<"INFO CHECK (one child image, get maps for the empty branches): ("<<node->id<<","<<dad->id<<")"<<endl;
            for(i=0; i<node_empty_branch_size; i++){
                
                ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part]=child_part_vec[0];
                ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = part_vec[0];
                
                if(back_branch_map){
                    child_part_vec[0]->link_neighbors.push_back(node->empty_br_node_nei[i]);
                    part_vec[0]->link_neighbors.push_back(node->empty_br_dad_nei[i]);
                    //cout<<i<<":"<<"("<<node->empty_br_node_nei[i]->node->id<<","<<node->empty_br_dad_nei[i]->node->id<<") -> ("<<child_part_vec[0]->node->id<<","<<part_vec[0]->node->id<<")"<<endl;
                }
            }
            node->empty_branches.clear();
            node->empty_br_node_nei.clear();
            node->empty_br_dad_nei.clear();
            
            if(back_taxon_map){
                //cout<<"Passing empty taxa to partition branches:"<<endl;
                for(i=0; i<node_empty_taxa_size; i++){
                    child_part_vec[0]->taxa_to_insert.push_back(node->empty_taxa[i]);
                    part_vec[0]->taxa_to_insert.push_back(node->empty_taxa[i]);
                }
                node->empty_taxa.clear();
            }
        }
        
        return;
    }
    
    // what's this below? -> When a subtree on the other side is empty. Do nothing for it.
    // In your case, you have have to map empty branches to (child_part_vec[0] ---- part_vec[1]) or (child_part_vec[1] ---- part_vec[0]). Does it matter, which nodes are mapped where? I suspect, that no. Aber there is always aber
    if (part_vec[0] == child_part_vec[1]) {

        //cout<<"CASE 2 SPECIAL: TWO CHILDREN HAVE IMAGEs, BUT THEY ARE THE SAME"<<endl<<endl;
        // ping-pong, out of sub-tree
        ASSERT(part_vec[1] == child_part_vec[0]);
        
        // Get maps for empty branches
        // the question is if your dad already has the details about all empty branches/taxa below ????? do not touch yet dad's empty branches, they should be considered at a later stage.
        // I think, in this situation you only need to map current branch, because node cannot have empty branches and taxa from its two children
        // again, i'm sure if the direction of neighbors is important or not...
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = child_part_vec[1];
        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            child_part_vec[1]->link_neighbors.push_back(dad_nei);
        }
        
        return;
    }
    
    //cout<<"CASE 2: TWO CHILDREN HAVE IMAGEs AND THEY ARE DIFFERENT"<<endl<<endl;
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
    
    if(back_branch_map){
        ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.push_back(nei);
        ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors.push_back(dad_nei);
    }

}

void Terrace::printMapInfo(){
    
    cout<<"Mapping info"<<endl<<endl;
    NodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    int part = 0;
    drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    for (vector<TerraceTree*>::iterator it = induced_trees.begin(); it != induced_trees.end(); it++, part++) {
        cout << endl << "Subtree for partition " << part << endl;
        // INFO: drawing of two-taxon tree fails.
        if((*it)->leafNum>2){
            (*it)->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE );
        }
        for (int i = 0; i < nodes1.size(); i++) {
            if(((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors.size()>0){
                Neighbor *nei1 = ((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part];
                Neighbor *nei2 = ((TerraceNeighbor*)nodes2[i]->findNeighbor(nodes1[i]))->link_neighbors[part];
                cout << nodes1[i]->findNeighbor(nodes2[i])->id << ":";
                if (nodes1[i]->isLeaf()) cout << nodes1[i]->name; else cout << nodes1[i]->id;
                cout << ",";
                if (nodes2[i]->isLeaf()) cout << nodes2[i]->name; else cout << nodes2[i]->id;
                //cout <<"("<<nodes1[i]->findNeighbor(nodes2[i])->length<<")"<< " -> ";
                cout<<" -> ";
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
                    //cout <<"("<<nei1->length<<")";
                }
                else cout << -1;
                cout << endl;
            }
        }
    }
    cout << endl;
}


void Terrace::printBackMapInfo(){

    // INFO: Also for upper level induced partition trees their induced subtrees are the same induced partition trees as for the master tree
    cout<<endl<<"Taxon mapping information:"<<endl;
    
    int i,j,k;
    NodeVector node_1, node_2;
    TerraceNeighbor *nei12, *nei21;
    
    for(i=0; i<part_num; i++){
        cout<<endl<<"Partition "<<i<<":"<<endl<<endl;
        node_1.clear();
        node_2.clear();
        induced_trees[i]->getBranches(node_1, node_2);
        for(j=0; j<node_1.size();j++){
            nei12 = (TerraceNeighbor*) node_1[j]->findNeighbor(node_2[j]);
            nei21 = (TerraceNeighbor*) node_2[j]->findNeighbor(node_1[j]);
            cout<<endl<<"* branch "<<nei12->id<<": "<<node_1[j]->id<<","<<node_2[j]->id<<endl;
            cout<<"+ link_neighbors:"<<endl;
            for(k=0; k<nei12->link_neighbors.size();k++){
                cout<<" - "<<k<<":"<<nei21->link_neighbors[k]->node->id<<","<<nei12->link_neighbors[k]->node->id<<endl;
            }
            cout<<endl<<"+ taxa:"<<endl;
            for(k=0; k<nei12->taxa_to_insert.size();k++){
                cout<<" - "<<k<<":"<<nei12->taxa_to_insert[k]->name<<" ("<<nei12->taxa_to_insert[k]->id<<")"<<endl;
            }
        }
    }
}

void Terrace::create_Top_Low_Part_Tree_Pairs(vector<Terrace*> &part_tree_pairs, Terrace *terrace){
    
    int i=0, j=0;
    NodeVector aux_taxon_nodes;
    vector<TerraceTree*> aux_induced_part_trees;
    IntVector parts;
    bool back_branch_map = false, back_taxon_map = true;
    
    for(i=0; i<terrace->part_num; i++){
        parts.clear();
        parts.push_back(i);
        aux_taxon_nodes.clear();
        
        // 1. get submatrix for taxa, which are common for the initial tree and for the top level induced partition tree. common taxa, are the leaves of the low level induced tree (assert: they should be with 1's!!!)
        induced_trees[i]->getTaxa(aux_taxon_nodes);
        PresenceAbsenceMatrix *aux_submatrix = new PresenceAbsenceMatrix();
        matrix->getSubPrAbMatrix(aux_taxon_nodes, aux_submatrix,&parts);
        
        // 2. now update matrix with taxa, which occur only on the top level induced partition tree and not on the initial tree. Simply add 0, because it does not occur on the common subtree, i.e. low level induced partition tree
        aux_taxon_nodes.clear();
        terrace->induced_trees[i]->getTaxa(aux_taxon_nodes);
        
        parts.clear(); // auxiliary use of the variable, just need an integer vec with 0 value.
        parts.push_back(0);
        
        for(j=0; j<aux_taxon_nodes.size(); j++){
            if(aux_submatrix->findTaxonID(aux_taxon_nodes[j]->name) == -1){
                aux_submatrix->pr_ab_matrix.push_back(parts);
                aux_submatrix->taxa_names.push_back(aux_taxon_nodes[j]->name);
                aux_submatrix->taxa_num +=1;
            }
        }
        
        //aux_submatrix->print_pr_ab_matrix();
        //terrace->induced_trees[i]->printTree(cout,WT_BR_LEN_ROUNDING + WT_NEWLINE);
        //cout<<endl;
        
        aux_induced_part_trees.clear();
        aux_induced_part_trees.push_back(induced_trees[i]);
        Terrace *aux_terrace = new Terrace(*(terrace->induced_trees[i]), aux_submatrix, aux_induced_part_trees);
        aux_terrace->linkTrees(back_branch_map, back_taxon_map);
        //aux_terrace->printMapInfo();
        //aux_terrace->printBackMapInfo();
        
        part_tree_pairs.push_back(aux_terrace);
        
    }
    
}
void Terrace::getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, vector<TerraceNeighbor*> *nei1_vec, vector<TerraceNeighbor*> *nei2_vec){
    
    int i, j, h;
    TerraceNode *node;
    TerraceNeighbor *nei, *dad_nei, *link_nei, *link_dad_nei;
    
    IntVector branch_ids, intersection;
    bool found = false;
    
    for(i=0;i<aux_terrace.size();i++){
        node = (TerraceNode*) aux_terrace[i]->findLeafName(taxon_name);
        if(node){
            assert(node->isLeaf());
            
            nei = (TerraceNeighbor*) node->neighbors[0];
            dad_nei = (TerraceNeighbor*) nei->node->findNeighbor(node);
            
            link_nei = (TerraceNeighbor*) nei->link_neighbors[0];
            link_dad_nei = (TerraceNeighbor*) dad_nei->link_neighbors[0];
            
            assert(link_nei->link_neighbors.size()==link_dad_nei->link_neighbors.size());
            
            if(branch_ids.empty()){
                for(j=0; j<link_nei->link_neighbors.size(); j++){
                    branch_ids.push_back(link_nei->link_neighbors[j]->id);
                    
                    // nei1_vec->push_back((TerraceNeighbor*)link_nei->link_neighbors[j]);
                    // nei2_vec->push_back((TerraceNeighbor*)link_dad_nei->link_neighbors[j]);
                }
            } else {
                
                // TODO: The way to keep track of allowed branches is not optimal, too many unnessary savings and clearings. Optimize at the next step. Also the intersection is O(n^2).
                nei1_vec->clear();
                nei2_vec->clear();
                
                // check, if the element occurs in branch ids
                for(j=0; j<link_nei->link_neighbors.size(); j++){
                    found = false;
                    for(h=0; h<branch_ids.size(); h++){
                        if(link_nei->link_neighbors[j]->id==branch_ids[h]){
                            found = true;
                            break;
                        }
                    }
                    if(found){
                        intersection.push_back(link_nei->link_neighbors[j]->id);
                        nei1_vec->push_back((TerraceNeighbor*)link_nei->link_neighbors[j]);
                        nei2_vec->push_back((TerraceNeighbor*)link_dad_nei->link_neighbors[j]);
                    }
                }
                
                branch_ids.clear();
                
                if(intersection.size()>0){
                    for(j=0; j<intersection.size(); j++){
                        branch_ids.push_back(intersection[j]);
                    }
                    intersection.clear();
                } else {
                    nei1_vec->clear();
                    nei2_vec->clear();
                    
                    // QUESTION: is this "break" breaks the FORLOOP with aux_terraces? in principle, if at some point intersection is empty, it means, that there is no point to look at other aux_terraces for allowed branches
                    break;
                }
            }
        }
    }
    
    /*if(!branch_ids.empty()){
        // get neighbours for branch with id branch_ids[i]
        for(i=0; i<branch_ids.size(); i++){
            
        }
    }*/
    
    
}

void Terrace::insertNewTaxon(string node_name, TerraceNeighbor *nei_1, TerraceNeighbor *nei_2){
    
    // TODO:
    // DONE: Node_1: get a corresponding node (it should be new one, create a new node)
    // DONE: Node_2: create a neighbour of the node
    // DONE: Update neighbours on the branch you need to insert the taxon to node_2 becomes a neighbour of noda_A and node_B, where the taxon is allowed to be inserted
    //       -> ABER: you have to provide a branch to be inserted to;
    //       -> QUESTION: how about updates of the link_neighbours? should you first update them and then update them on the final tree
    // Update mapping from parent tree to low-level induced partition trees
    // Update mapping from top-level to low-level trees (only if the top-level has the corresponding taxon)
    // repeat: get allowed positions for the next taxon
   
    // IMPORTANT!!!!!! WARNING: CAREFULL about calling a presence_absence matrix entries by the leaf id!!!!!! CHECK the corresponding code, because your maps will be fucked up!
    TerraceNode *node_1, *node_2;
    
    newNode(nodeNum, node_name.c_str());
    node_1 = (TerraceNode*)(this->findNodeID(nodeNum));
    
    leafNum += 1;
    nodeNum += 1;
    taxa_num += 1;

    newNode(nodeNum);
    node_2 = (TerraceNode*)(this->findNodeID(nodeNum));
    nodeNum += 1;
    branchNum += 1;
    
    node_1->addNeighbor(node_2, -1);
    
    nei_1->node->updateNeighbor(nei_2->node, node_2, -1);
    nei_2->node->updateNeighbor(nei_1->node, node_2, -1);
    
    node_2->addNeighbor(node_1, -1);
    node_2->addNeighbor(nei_1->node, -1);
    node_2->addNeighbor(nei_2->node, -1);
    
}

// INFO: TWO FUNCTIONs BELOW ARE not finished (kind of obvious), but also maybe unnecessary in the current setting of the mapping
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

void Terrace::clearEmptyBranchAndTaxaINFO(TerraceNode *node, TerraceNode *dad){
    
    if (!node) {
        if (!root->isLeaf())
            node = (TerraceNode*) root;
        else
            node = (TerraceNode*) root->neighbors[0]->node;
        ASSERT(node);
    }
    
    if(node->empty_taxa.size()>0){
        node->empty_taxa.clear();
    }
    
    if(node->empty_branches.size()>0){
        node->empty_branches.clear();
        node->empty_br_node_nei.clear();
        node->empty_br_dad_nei.clear();
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<endl;
        clearEmptyBranchAndTaxaINFO((TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
}

void Terrace::cleanAllLinkNeighboursAndTaxa(){
 
    cleanAllLinkINFO();
    
}
