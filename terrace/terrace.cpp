//
//  terrace.cpp
//  iqtree
//
//  Created by Olga on 14.09.20.
//

#include "terrace.hpp"
#include "terracenode.hpp"
#include "tree/mtreeset.h"
#include "utils/timeutil.h"

Terrace::Terrace(){};
Terrace::~Terrace(){
    for(auto it=induced_trees.begin(); it<induced_trees.end();it++){
        delete (*it);
    }
    induced_trees.clear();
    
    delete matrix;
};

void Terrace::init(){
    
    taxa_num = 0;
    part_num = 0;
    intermediated_trees_num = 0;
    terrace_trees_num = 0;
    dead_ends_num = 0;
    terrace_out = true;
    trees_out_lim = 0;
    
    terrace_max_trees = 1000000;
    intermediate_max_trees = 10000000;
    seconds_max = 604800; // the default value is 7 days
    
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
    
    //renameTaxa();
    
    get_part_trees();
    //printInfo();
};

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix *m){
    
    init();

    if(tree.leafNum>2){
        copyTree(&tree);
    }else{
        string taxa_set = "";
        vector<uint32_t> check_int;
        check_int.resize(2,1);
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        copyTree(&tree,taxa_set);
    }
    
    matrix = m;
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();
    
    NodeVector taxa_nodes;
    getTaxa(taxa_nodes);
    //printTree(cout);
    
    matrix->reorderAccordingToTree(taxa_nodes);
    //matrix->print_pr_ab_matrix();
    
    get_part_trees();
    //printInfo();
}

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix *m, vector<TerraceTree*> input_induced_trees){
 
    init();
    
    if(tree.leafNum>2){
        copyTree(&tree);
    }else{
        string taxa_set = "";
        vector<uint32_t> check_int;
        check_int.resize(2,1);
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        copyTree(&tree,taxa_set);
    }
    matrix = m;
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();
    
    NodeVector taxa_nodes;
    getTaxa(taxa_nodes);
    //printTree(cout);
    
    matrix->reorderAccordingToTree(taxa_nodes);
    //matrix->print_pr_ab_matrix();
    
    set_part_trees(input_induced_trees);
    //printInfo();
    
}

void Terrace::get_part_trees(){
    
    assert(matrix!=nullptr && "ERROR: Presence-absence matrix is absent! I cannot get induced partition trees..");
    //cout<<"Preparing induced partition trees..."<<endl;
    
    int id,k;
    string taxa_set = "";
    string taxon_name;
    NodeVector taxa_nodes;
    
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(taxa_num);
    
    getTaxa(taxa_nodes);
    int leafNum_part_tree = 0;
    for(int part=0; part<part_num; part++){
        leafNum_part_tree = 0;
        TerraceTree* induced_tree = new TerraceTree();
        for(it2=taxa_nodes.begin();it2!=taxa_nodes.end();it2++){
            taxon_name=(*it2)->name;
            id=matrix->findTaxonID(taxon_name);
            assert(id != -1 && "Not all of the taxa appear in the pr_ab_matrix!");
            check_int[id] = matrix->pr_ab_matrix[id][part];
            //cout<<"Taxon["<<id<<"] in partition "<<part<<" is "<<check_int[id]<<endl;
            if(check_int[id]==1){
                leafNum_part_tree+=1;
            }
        }
        taxa_set.clear();
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        
        if(leafNum_part_tree>1){
            induced_tree->copyTree(this,taxa_set);
        } else if(leafNum_part_tree==1){
            for(k=0; k<taxa_num; k++){
                if(check_int[k]==1){
                    taxon_name = matrix->taxa_names[k];
                    break;
                }
            }
            
            induced_tree->root = induced_tree->newNode(0, taxon_name.c_str());
            induced_tree->leafNum = 1;
            induced_tree->nodeNum = 1;
        }
        
        induced_trees.push_back(induced_tree);
    }
}

void Terrace::set_part_trees(vector<TerraceTree*> input_induced_trees){
    
    for(int i=0; i<input_induced_trees.size(); i++){
        induced_trees.push_back(input_induced_trees[i]);
    }
    
    // TODO: you need a check, that your input_induced_trees indeed correspond to presence-absence pattern of pr_ab_matrix
}

void Terrace::printInfo(){

    cout<<endl<<"=================================================="<<endl<<"Printing terrace information:"<<endl<<endl;
    cout<<"Terrace representative tree:"<<endl;
    print_terrace_tree();
    cout<<endl;
    
    matrix->print_pr_ab_matrix();
    
    if(induced_trees.size()>0){
        int i=0;
        cout<<"Induced partition trees:"<<endl;
        for(vector<TerraceTree*>::iterator it = induced_trees.begin(); it < induced_trees.end(); it++){
            i++;
            cout<<"Part["<<i<<"]: ";
            (*it)->print_terrace_tree(false);
        }
    }
    cout<<endl<<"=================================================="<<endl<<endl;
}

void Terrace::linkTrees(bool back_branch_map, bool back_taxon_map){
    NodeVector part_taxa;
    
    for(int part=0; part<part_num; part++){
        if(induced_trees[part]->leafNum>2){
            part_taxa.clear();
            matrix->getPartTaxa(part, this, induced_trees[part], part_taxa);
            //clearEmptyBranchAndTaxaINFO();
            linkTree(part, part_taxa, back_branch_map, back_taxon_map);
        }
    }
    
}

void Terrace::linkTree(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad){
    
    // SEHR WICHTIG! WARNING: do not mix mapping from the parent tree and upper level induced partition trees, because empty branches and empty taxa will be messed up on partition trees.
    
    // INFO/CHECK: make sure that your branches have proper unique ids, further code relies on that.
    
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
        // DONE: FIXED. BUG: I think this is the case. again issue with the node ids. Find another way of getting node->id. Check if this is actually what causes problem
        //TerraceNode *node_part = (TerraceNode*)part_taxa[node->id]; //part_taxa[] is a node on PARTITION tree
        TerraceNode *node_part = nullptr;
        for(NodeVector::iterator it=part_taxa.begin(); it<part_taxa.end(); it++){
            if((*it)){
                if(node->name == (*it)->name){
                    node_part = (TerraceNode*)(*it);
                    break;
                }
            }
        }
        
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
            
            if(back_taxon_map){
                // Back map from induced partition tree onto parent tree. Here used for low to top induced partition trees backward map
                // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors_lowtop_back.push_back(nei);
                dad_part_nei->link_neighbors_lowtop_back.push_back(dad_nei);
            }
            
        } else {
            
            // the empty branches are needed for a map from parent to partition
            //dad->empty_branches.push_back(nei->id);
            dad->empty_br_dad_nei.push_back(dad_nei);
            dad->empty_br_node_nei.push_back(nei);
            
            // INFO: The actual taxon maps are not needed, but they are convenient for the check, but at the current state with local updates, will be incorrect
            /*if(back_taxon_map){
                // the empty taxa are needed for a map from induced partition tree to "common" induced partition tree
                // when you know the link_neighbor, you have to map all empty taxa to that neighbor
                dad->empty_taxa.push_back(node);
                //cout<<"INFO CHECK: the length of the empty leaves: dad("<<dad->id<<") = "<<dad->empty_taxa.size()<<endl;
            }*/
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
        //if(node->empty_taxa.size()>0 or node->empty_branches.size()>0){
        if(!node->empty_br_node_nei.empty()){
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
                        
                        if(back_taxon_map){
                            node_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_dad_nei[i]);
                        }
                        
                    }
                    
                    /*if(back_taxon_map){
                        //cout<<endl<<"INFO CHECK (in IF!dad): node="<<node->id<<" | empty taxa:"<<endl;
                        for(i=0; i<node->empty_taxa.size(); i++){
                            //cout<<" "<<node->empty_taxa[i]->id<<",";
                            node_nei_part->taxa_to_insert.push_back(node->empty_taxa[i]);
                            dad_nei_part->taxa_to_insert.push_back(node->empty_taxa[i]);
                        }
                        //cout<<endl;
                        node->empty_taxa.clear();
                    }*/
                        
                    //node->empty_branches.clear();
                    node->empty_br_node_nei.clear();
                    node->empty_br_dad_nei.clear();
                    
                    return;
                }
            }
            
            // WARNING: I think, it should never rich this point, because at least one neighbour should have a corresponding branch on partition tree.
            cout<<"ERROR: hey, the assumption was, that the code should never come here..."<<endl;
        }
        return;
    }
    linkBranch(part, nei, dad_nei, back_branch_map, back_taxon_map);
}

void Terrace::linkBranch(int part, TerraceNeighbor *nei, TerraceNeighbor *dad_nei, bool back_branch_map, bool back_taxon_map) {
    
    TerraceNode *node = (TerraceNode*)dad_nei->node;
    TerraceNode *dad = (TerraceNode*)nei->node;
    nei->link_neighbors[part] = NULL;
    dad_nei->link_neighbors[part] = NULL;
    vector<TerraceNeighbor*> part_vec;
    vector<TerraceNeighbor*> child_part_vec;
    
    //cout<<"Linking branch ("<<node->id<<","<<dad->id<<");"<<endl;
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        if (((TerraceNeighbor*)*it)->link_neighbors[part]) {
            part_vec.push_back((TerraceNeighbor*)(((TerraceNeighbor*)*it)->link_neighbors[part]));
            child_part_vec.push_back((TerraceNeighbor*)(((TerraceNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]));
            
            // INFO: "folding" of branches, which can be folded to branch and map to "nodeA-nodeB" or "nodeB-nodeA", children neighbors can be different. Therefore, in some cases the assertion won't be fulfilled, but this is not an error and, thus, the assertion is commented out.
            //ASSERT(child_part_vec.back()->node == child_part_vec.front()->node || child_part_vec.back()->id == child_part_vec.front()->id);
        }
    }
    
    if(child_part_vec.size()==2){
        if(child_part_vec[0]->node != child_part_vec[1]->node){
            //cout<<"My weird case I..."<<endl;
            //cout<<"child_part_vec[0]="<<child_part_vec[0]->node->id<<endl;
            //cout<<"child_part_vec[1]="<<child_part_vec[1]->node->id<<endl;
            //cout<<"part_vec[0]="<<part_vec[0]->node->id<<endl;
            //cout<<"part_vec[1]="<<part_vec[1]->node->id<<endl;
            TerraceNeighbor * nei_aux;
            if(part_vec[0]->node == part_vec[1]->node){
                nei_aux = part_vec[0];
                part_vec[0]=child_part_vec[0];
                child_part_vec[0]=nei_aux;
                nei_aux=part_vec[1];
                part_vec[1]=child_part_vec[1];
                child_part_vec[1]=nei_aux;
            }else{
                for(int i=0; i<part_vec.size(); i++){
                    for(int j=0; j<child_part_vec.size(); j++){
                        if(part_vec[i]->node == child_part_vec[j]->node and part_vec[1-i]->node != child_part_vec[1-j]->node){
                            nei_aux = part_vec[i];
                            part_vec[i] = child_part_vec[1-j];
                            child_part_vec[1-j] = nei_aux;
                            break;
                        }
                    }
                }
            }
            //cout<<"child_part_vec[0]="<<child_part_vec[0]->node->id<<endl;
            //cout<<"child_part_vec[1]="<<child_part_vec[1]->node->id<<endl;
            //cout<<"part_vec[0]="<<part_vec[0]->node->id<<endl;
            //cout<<"part_vec[1]="<<part_vec[1]->node->id<<endl;
        }
    }

    
    
    // QUESTION: node or dad empty branch/taxa?
    // QUESTION: if you start mapping the dad, you might encounter situations, when a dad does not know yet, that there are some empty branches.
    // So at the moment you need only empty branches of node and the dad will be mapped at the final state (i.e. final node - the start of the tree traversal)
    
    //int node_empty_branch_size = node->empty_branches.size();
    //int node_empty_taxa_size = node->empty_taxa.size();
    
    if (part_vec.empty()){
        int i=0;
        //cout<<"CASE 0: CHILDREN DO NOT HAVE IMAGE"<<endl<<endl;
        if(!node->empty_br_dad_nei.empty()){
            for(i=0; i<node->empty_br_dad_nei.size(); i++){
                //dad->empty_branches.push_back(node->empty_branches[i]);
                dad->empty_br_dad_nei.push_back(node->empty_br_dad_nei[i]);
                dad->empty_br_node_nei.push_back(node->empty_br_node_nei[i]);
            }
            //node->empty_branches.clear();
            node->empty_br_node_nei.clear();
            node->empty_br_dad_nei.clear();
        }
        
        // since below node the subtrees are empty, add also current branch to an empty list of the dad
        //dad->empty_branches.push_back(nei->id);
        dad->empty_br_dad_nei.push_back(dad_nei);
        dad->empty_br_node_nei.push_back(nei);
        
        /*if(back_taxon_map){
            if(node_empty_taxa_size>0){
                cout<<"INFO CHECK (in two children nodes are empty): passing empty taxa to dad="<<dad->id<<" :"<<endl;
                if(dad->empty_taxa.size()>0){
                    cout<<"currently list of empty taxa for dad:";
                    for(i=0; i<dad->empty_taxa.size(); i++){
                        cout<<" "<<dad->empty_taxa[i]->id;
                    }
                    cout<<endl;
                }
         
                
                for(i=0; i<node_empty_taxa_size; i++){
                    //cout<<" "<<node->empty_taxa[i]->id;
                    dad->empty_taxa.push_back(node->empty_taxa[i]);
                }
                //cout<<endl;
                node->empty_taxa.clear();
            }
        }*/
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
        
        if(back_taxon_map){
            child_part_vec[0]->link_neighbors_lowtop_back.push_back(nei);
            part_vec[0]->link_neighbors_lowtop_back.push_back(dad_nei);
        }
        
        // Get maps for empty branches
        if(node->empty_br_node_nei.size()>0){
            //cout<<"INFO CHECK (one child image, get maps for the empty branches): ("<<node->id<<","<<dad->id<<")"<<endl;
            for(i=0; i<node->empty_br_node_nei.size(); i++){
                
                ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part]=child_part_vec[0];
                ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = part_vec[0];
                
                if(back_branch_map){
                    child_part_vec[0]->link_neighbors.push_back(node->empty_br_node_nei[i]);
                    part_vec[0]->link_neighbors.push_back(node->empty_br_dad_nei[i]);
                    //cout<<i<<":"<<"("<<node->empty_br_node_nei[i]->node->id<<","<<node->empty_br_dad_nei[i]->node->id<<") -> ("<<child_part_vec[0]->node->id<<","<<part_vec[0]->node->id<<")"<<endl;
                }
                
                if(back_taxon_map){
                    child_part_vec[0]->link_neighbors_lowtop_back.push_back(node->empty_br_node_nei[i]);
                    part_vec[0]->link_neighbors_lowtop_back.push_back(node->empty_br_dad_nei[i]);
                }
                
            }
            //node->empty_branches.clear();
            node->empty_br_node_nei.clear();
            node->empty_br_dad_nei.clear();
            
            /*if(back_taxon_map){
                //cout<<"Passing empty taxa to partition branches:"<<endl;
                for(i=0; i<node_empty_taxa_size; i++){
                    child_part_vec[0]->taxa_to_insert.push_back(node->empty_taxa[i]);
                    part_vec[0]->taxa_to_insert.push_back(node->empty_taxa[i]);
                }
                node->empty_taxa.clear();
            }*/
        }
        return;
    }
    
    // what's this below? -> When a subtree on the other side is empty. Do nothing for it.
    // In your case, you have have to map empty branches to (child_part_vec[0] ---- part_vec[1]) or (child_part_vec[1] ---- part_vec[0]). Does it matter, which nodes are mapped where? I suspect, that no. Aber there is always aber
    if (part_vec[0] == child_part_vec[1]) {
        // ping-pong, out of sub-tree
        ASSERT(part_vec[1] == child_part_vec[0]);
        
        // Get maps for empty branches
        // the question is if your dad already has the details about all empty branches/taxa below ????? do not touch yet dad's empty branches, they should be considered at a later stage.
        // I think, in this situation you only need to map current branch, because node cannot have empty branches and taxa from its two children
        // again, i'm not sure if the direction of neighbors is important or not...
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = child_part_vec[1];

        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            child_part_vec[1]->link_neighbors.push_back(dad_nei);
        }
        if(back_taxon_map){
            child_part_vec[0]->link_neighbors_lowtop_back.push_back(nei);
            child_part_vec[1]->link_neighbors_lowtop_back.push_back(dad_nei);
        }
        return;
    }
    
    // WARNING/INFO: FOR THE UPDATE_MAP: CASE 2 (two identical images for child branches), might also happen, when
    // child_part_vec[0]=child_part_vec[1]
    // part_vec[0]=part_vec[1]
    // this is due to the fact, that all branches have a link_neighbors and the direction in which the branch "folds" is not certain
    //cout<<"DEBUGGING:"<<endl;
    //cout<<"child_part_vec[0]="<<child_part_vec[0]->node->id<<endl;
    //cout<<"child_part_vec[1]="<<child_part_vec[1]->node->id<<endl;
    //cout<<"part_vec[0]="<<part_vec[0]->node->id<<endl;
    //cout<<"part_vec[1]="<<part_vec[1]->node->id<<endl;
    
    if(child_part_vec[0]==child_part_vec[1]){
        //cout<<"My weird case II..."<<endl;
        
        assert(part_vec[0]==part_vec[1]);
        
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = part_vec[0];
        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            part_vec[0]->link_neighbors.push_back(dad_nei);
        }
        if(back_taxon_map){
            child_part_vec[0]->link_neighbors_lowtop_back.push_back(nei);
            part_vec[0]->link_neighbors_lowtop_back.push_back(dad_nei);
        }
        return;
    }
    
    
    //cout<<"CASE 2: TWO CHILDREN HAVE IMAGEs AND THEY ARE DIFFERENT"<<endl<<endl;
    TerraceNode *node_part = (TerraceNode*) child_part_vec[0]->node;
    TerraceNode *dad_part = nullptr;
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
    
    if(back_taxon_map){
        ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors_lowtop_back.push_back(nei);
        ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors_lowtop_back.push_back(dad_nei);
    }
}

void Terrace::update_map(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad){
    
    //if(!dad){
    //    cout<<"++++++++++++ PARTITION "<<part<<"++++++++++++"<<endl;
    //    printMapInfo(3);
    //}
    
    //if(part==0){
    //    printMapInfo();
    //}
    
    /*
     For insertion:
         1. You need some kind of preprocessing step to remove link_neighbours of branches, which were allowed for the inserted taxon.
         2. Then call this function recursively.
         3. What shall we do with a drunken sailor early in the morning?:) take into account cases when the inserted taxon is 1st or 2nd for some partition.
     */
    
    if (node->isLeaf()){
        if(!dad){
            dad = (TerraceNode*)node->neighbors[0]->node;
            if(!dad->isLeaf()){
                node = dad;
                dad = nullptr;
            }
        }
    }
    
    if(!dad){
        FOR_NEIGHBOR_DECLARE(node, dad, it) {
            //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<endl;
            update_map(part, part_taxa, back_branch_map, back_taxon_map, (TerraceNode*) (*it)->node, (TerraceNode*) node);
        }
        
        // Check, if the final dad has empty_branches and empty_taxa and map them to branch, which is available.
        // Note, that if there are some empty branches/taxa, there will be exactly one branch available for mapping.
        //if(node->empty_taxa.size()>0 or node->empty_branches.size()>0){
        if(!node->empty_br_node_nei.empty()){
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
                        
                        if(back_taxon_map){
                            node_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_dad_nei[i]);
                        }
                    }
                    
                    /*if(back_taxon_map){
                        //cout<<endl<<"INFO CHECK (in IF!dad): node="<<node->id<<" | empty taxa:"<<endl;
                        for(i=0; i<node->empty_taxa.size(); i++){
                            //cout<<" "<<node->empty_taxa[i]->id<<",";
                            node_nei_part->taxa_to_insert.push_back(node->empty_taxa[i]);
                            dad_nei_part->taxa_to_insert.push_back(node->empty_taxa[i]);
                        }
                        //cout<<endl;
                        node->empty_taxa.clear();
                    }*/
                    
                    //node->empty_branches.clear();
                    node->empty_br_node_nei.clear();
                    node->empty_br_dad_nei.clear();
                    
                    return;
                }
            }
        }
        return;
        
    } else {
        TerraceNeighbor *nei = (TerraceNeighbor*)node->findNeighbor(dad);
        TerraceNeighbor *dad_nei = (TerraceNeighbor*)dad->findNeighbor(node);
    
        // if there is a map already, return
        if(nei->link_neighbors[part]){
            if(dad_nei->link_neighbors[part])
                return;
        } else {
    
            if (node->isLeaf()) {
                ASSERT(dad);
                TerraceNode *node_part = nullptr;
            
                for(NodeVector::iterator it=part_taxa.begin(); it<part_taxa.end(); it++){
                    if((*it)){
                        if(node->name == (*it)->name){
                            node_part = (TerraceNode*)(*it);
                            break;
                        }
                    }
                }
            
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
                    
                    if(back_taxon_map){
                        // Back map from induced partition tree onto parent tree (for LOW-TOP backmap)
                        // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                        ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors_lowtop_back.push_back(nei);
                        dad_part_nei->link_neighbors_lowtop_back.push_back(dad_nei);
                    }
                    
                } else {
                    // the empty branches are needed for a map from parent to partition
                    //dad->empty_branches.push_back(nei->id);
                    dad->empty_br_dad_nei.push_back(dad_nei);
                    dad->empty_br_node_nei.push_back(nei);
                    
                    /*if(back_taxon_map){
                        // the empty taxa are needed for a map from induced partition tree to "common" induced partition tree
                        // when you know the link_neighbor, you have to map all empty taxa to that neighbor
                        dad->empty_taxa.push_back(node);
                        //cout<<"INFO CHECK: the length of the empty leaves: dad("<<dad->id<<") = "<<dad->empty_taxa.size()<<endl;
                    }*/
                }
                return;
            }
        
            //cout<<"MASTER NODE:"<<node->id<<endl;
            FOR_NEIGHBOR_DECLARE(node, dad, it) {
                //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<endl;
                update_map(part, part_taxa, back_branch_map, back_taxon_map, (TerraceNode*) (*it)->node, (TerraceNode*) node);
            }

            linkBranch(part, nei, dad_nei, back_branch_map, back_taxon_map);
        }
    }
    
}

void Terrace::printMapInfo(int partition){
    
    cout<<"Mapping info"<<endl<<endl;
    NodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    int part = 0;
    drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    if(partition==-1){
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
    } else {
        part = partition;
        cout << endl << "Subtree for partition " << part << endl;
        // INFO: drawing of two-taxon tree fails.
        if(induced_trees[part]->leafNum>2){
            induced_trees[part]->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE );
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
    cout<<endl<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<endl;
    cout<<endl<<"BACKWARD mapping information:"<<endl;
    cout<<endl<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<endl;
    
    int i,j,k;
    NodeVector node_1, node_2;
    TerraceNeighbor *nei12, *nei21;
    
    i=71;
    //for(i=0; i<part_num; i++){
        cout<<endl<<"---------------------------------------------"<<endl;
        cout<<endl<<"Partition "<<i<<":"<<endl;
        cout<<endl<<"---------------------------------------------"<<endl;
        if(induced_trees[i]->leafNum>1){
            node_1.clear();
            node_2.clear();
            induced_trees[i]->getBranches(node_1, node_2);
            for(j=0; j<node_1.size();j++){
                nei12 = (TerraceNeighbor*) node_1[j]->findNeighbor(node_2[j]);
                nei21 = (TerraceNeighbor*) node_2[j]->findNeighbor(node_1[j]);
                cout<<endl<<"* branch "<<nei12->id<<": ";
                if(node_1[j]->isLeaf()){
                    cout<<node_1[j]->name;
                }
                cout<<"("<<node_1[j]->id<<")"<<",";
                if(node_2[j]->isLeaf()){
                    cout<<node_2[j]->name;
                }
                cout<<"("<<node_2[j]->id<<")"<<endl;
                cout<<"+ link_neighbors:"<<endl;
                if(!nei12->link_neighbors.empty()){
                    for(k=0; k<nei12->link_neighbors.size();k++){
                        cout<<" - "<<nei21->link_neighbors[k]->id<<":";
                        if(nei21->link_neighbors[k]->node->isLeaf()){
                            cout<<nei21->link_neighbors[k]->node->name;
                        }
                        cout<<"("<<nei21->link_neighbors[k]->node->id<<"),";
                        if(nei12->link_neighbors[k]->node->isLeaf()){
                            cout<<nei12->link_neighbors[k]->node->name;
                        }
                        cout<<"("<<nei12->link_neighbors[k]->node->id<<")"<<endl;
                    }
                }
                //cout<<endl<<"+ taxa:"<<endl;
                //for(k=0; k<nei12->taxa_to_insert.size();k++){
                //    cout<<" - "<<k<<":"<<nei12->taxa_to_insert[k]->name<<" ("<<nei12->taxa_to_insert[k]->id<<")"<<endl;
                //}
                cout<<"+ link_neighbors_lowtop_back:"<<endl;
                if(!nei12->link_neighbors_lowtop_back.empty()){
                    for(k=0; k<nei12->link_neighbors_lowtop_back.size();k++){
                        cout<<" - "<<nei21->link_neighbors_lowtop_back[k]->id<<":";
                        if(nei21->link_neighbors_lowtop_back[k]->node->isLeaf()){
                            cout<<nei21->link_neighbors_lowtop_back[k]->node->name;
                        }
                        cout<<"("<<nei21->link_neighbors_lowtop_back[k]->node->id<<"),";
                        if(nei12->link_neighbors_lowtop_back[k]->node->isLeaf()){
                            cout<<nei12->link_neighbors_lowtop_back[k]->node->name;
                        }
                        cout<<"("<<nei12->link_neighbors_lowtop_back[k]->node->id<<")"<<endl;
                    }
                }
            }
        }
   // }
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
        PresenceAbsenceMatrix *aux_submatrix = new PresenceAbsenceMatrix();
        if(induced_trees[i]->leafNum>0){
            induced_trees[i]->getTaxa(aux_taxon_nodes);
            matrix->getSubPrAbMatrix(aux_taxon_nodes, aux_submatrix,&parts);
        }
        // 2. now update matrix with taxa, which occur only on the top level induced partition tree and not on the initial tree. Simply add 0, because it does not occur on the common subtree, i.e. low level induced partition tree
        aux_taxon_nodes.clear();
        terrace->induced_trees[i]->getTaxa(aux_taxon_nodes);
        
        parts.clear(); // auxiliary use of the variable, just need an integer vec with 0 value.
        parts.push_back(0);
        
        for(j=0; j<aux_taxon_nodes.size(); j++){
            if(induced_trees[i]->leafNum>0){
                if(aux_submatrix->findTaxonID(aux_taxon_nodes[j]->name) == -1){
                    aux_submatrix->pr_ab_matrix.push_back(parts);
                    aux_submatrix->taxa_names.push_back(aux_taxon_nodes[j]->name);
                    aux_submatrix->taxa_num +=1;
                }
            }else{
                aux_submatrix->pr_ab_matrix.push_back(parts);
                aux_submatrix->taxa_names.push_back(aux_taxon_nodes[j]->name);
                aux_submatrix->taxa_num +=1;
            }
        }
        aux_submatrix->part_num = 1;
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
        
        //aux_terrace->matrix->print_pr_ab_matrix();
    }
}

void Terrace::getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, NodeVector *node1_vec, NodeVector *node2_vec){
    // INFO/CHECK: make sure that your branches have proper unique ids, below code relies on that.
    /*cout<<endl<<"**********************************************"<<endl<<endl;
    cout<<"IN getALLOWEDbranches"<<endl;
    cout<<endl<<"**********************************************"<<endl<<endl;

    printMapInfo();
    printBackMapInfo();*/
    
    int i, j, h;
    TerraceNode *node;
    TerraceNeighbor *nei, *dad_nei, *link_nei, *link_dad_nei;
    
    IntVector branch_ids, intersection;
    bool found = false;
    
    for(i=0;i<aux_terrace.size();i++){
        node = (TerraceNode*) aux_terrace[i]->findLeafName(taxon_name);
        if(node){
            assert(node->isLeaf());
            
            if(induced_trees[i]->leafNum>2){
                nei = (TerraceNeighbor*) node->neighbors[0];
                dad_nei = (TerraceNeighbor*) nei->node->findNeighbor(node);
                
                link_nei = (TerraceNeighbor*) nei->link_neighbors[0];
                link_dad_nei = (TerraceNeighbor*) dad_nei->link_neighbors[0];
                
                assert(link_nei->link_neighbors.size()==link_dad_nei->link_neighbors.size());
            }
            if(branch_ids.empty()){
                if(induced_trees[i]->leafNum>2){
                    for(j=0; j<link_nei->link_neighbors.size(); j++){
                        branch_ids.push_back(link_nei->link_neighbors[j]->id);
                        
                        intersection.push_back(link_nei->link_neighbors[j]->id);
                        node1_vec->push_back((TerraceNode*)link_nei->link_neighbors[j]->node);
                        node2_vec->push_back((TerraceNode*)link_dad_nei->link_neighbors[j]->node);
                    }
                } else {
                    //cout<<"All branches are allowed for this tree. Push them back"<<endl;
                    NodeVector branch_end_1, branch_end_2;
                    Neighbor* aux_nei;
                    this->getBranches(branch_end_1, branch_end_2);
                    for(int k=0; k<branch_end_1.size(); k++){
                        aux_nei = branch_end_1[k]->findNeighbor(branch_end_2[k]);
                        
                        branch_ids.push_back(aux_nei->id);
                        intersection.push_back(aux_nei->id);
                        
                        node1_vec->push_back((TerraceNode*)branch_end_1[k]);
                        node2_vec->push_back((TerraceNode*)branch_end_2[k]);
                    }
                }
            } else if(induced_trees[i]->leafNum>2){ // otherwise, there is no constraint
                
                // TODO: The way to keep track of allowed branches is not optimal, too many unnessary savings and clearings. Optimize in the next step. Also the intersection is O(n^2).
                //nei1_vec->clear();
                //nei2_vec->clear();
                intersection.clear();
                node1_vec->clear();
                node2_vec->clear();
                
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
                        //nei1_vec->push_back((TerraceNeighbor*)link_nei->link_neighbors[j]);
                        //nei2_vec->push_back((TerraceNeighbor*)link_dad_nei->link_neighbors[j]);
                        node1_vec->push_back((TerraceNode*)link_nei->link_neighbors[j]->node);
                        node2_vec->push_back((TerraceNode*)link_dad_nei->link_neighbors[j]->node);
                    }
                }
                
                branch_ids.clear();
                
                if(intersection.size()>0){
                    for(j=0; j<intersection.size(); j++){
                        branch_ids.push_back(intersection[j]);
                    }
                    intersection.clear();
                } else {
                    //nei1_vec->clear();
                    //nei2_vec->clear();
                    
                    node1_vec->clear();
                    node2_vec->clear();
                    
                    // QUESTION: is this "break" breaks the FORLOOP with aux_terraces? in principle, if at some point intersection is empty, it means, that there is no point to look at other aux_terraces for allowed branches
                    break;
                }
            }
        }
    }
}

void Terrace::extendNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> part_tree_pairs,vector<int> pr_ab_info){
    
    TerraceNeighbor *nei_1, *nei_2;
    nei_1 = (TerraceNeighbor*) node_1_branch->findNeighbor(node_2_branch);
    nei_2 = (TerraceNeighbor*) node_2_branch->findNeighbor(node_1_branch);
    
    NodeVector induced_part_tree_branch_1, induced_part_tree_branch_2;
    induced_part_tree_branch_1.resize(part_num,nullptr);
    induced_part_tree_branch_2.resize(part_num,nullptr);
    TerraceNeighbor *nei_part_1, *nei_part_2;
    NodeVector part_taxa;
    
    int i,j;

    for(i=0; i<part_num; i++){
        
        if(induced_trees[i]->leafNum>2){
            
            induced_part_tree_branch_1[i]=((TerraceNeighbor*)nei_2->link_neighbors[i])->node;
            induced_part_tree_branch_2[i]=((TerraceNeighbor*)nei_1->link_neighbors[i])->node;
            
            nei_part_1 = (TerraceNeighbor*)induced_part_tree_branch_1[i]->findNeighbor(induced_part_tree_branch_2[i]);
            nei_part_2 = (TerraceNeighbor*)induced_part_tree_branch_2[i]->findNeighbor(induced_part_tree_branch_1[i]);
        
            if(part_tree_pairs[i]->findLeafName(node_name)){
                
                // --------------- CLEARING MAPS FROM AGILE TREE TO COMMON PART SUBTREES and BACK --------------------------------------------------------------------------
                // STEP 1: clearing pointers for all branches that mapped to the involved branch on induced partition tree (clearing forward maps)
                for(j=0; j<nei_part_1->link_neighbors.size();j++){
                    ((TerraceNeighbor*)nei_part_1->link_neighbors[j])->link_neighbors[i] = nullptr;
                    ((TerraceNeighbor*)nei_part_2->link_neighbors[j])->link_neighbors[i] = nullptr;
                }
                
                // STEP 2: clearing backward maps from induced partition trees to parent tree
                nei_part_1->link_neighbors.clear();
                nei_part_2->link_neighbors.clear();
                
                
                // --------------- CLEARING MAPS FROM INDUCED PART TREE TO COMMON PART SUBTREES and BACK --------------------------------------------------------------------------
                // STEP 1: clearing pointers for all branches that mapped to the involved branch on induced partition tree (clearing forward maps)
                for(j=0; j<nei_part_1->link_neighbors_lowtop_back.size();j++){
                    ((TerraceNeighbor*)nei_part_1->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                    ((TerraceNeighbor*)nei_part_2->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                }
                
                // STEP 2: clearing backward maps from induced partition trees to parent tree
                nei_part_1->link_neighbors_lowtop_back.clear();
                nei_part_2->link_neighbors_lowtop_back.clear();
                
                //part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa();

                // --------------- INSERT A TAXON ON COMMON PARTITION SUBTREE ----------------------------------------------------------------------------------------------------
                induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) induced_part_tree_branch_1[i], (TerraceNode*) induced_part_tree_branch_2[i]);
                part_tree_pairs[i]->matrix->pr_ab_matrix[part_tree_pairs[i]->matrix->findTaxonID(node_name)][0]=1;
            

                // --------------- UPDATE MAPs LOCALLY for TOP LEVEL INDUCED PART TREES ------------------------------------------------------------------------------------------
                part_taxa.clear();
                part_tree_pairs[i]->matrix->getPartTaxa(0, part_tree_pairs[i], part_tree_pairs[i]->induced_trees[0], part_taxa);
                
                //part_tree_pairs[i]->linkTree(0, part_taxa, false, true);
                TerraceNode *central_node_part = (TerraceNode*) part_tree_pairs[i]->findLeafName(node_name)->neighbors[0]->node;
                part_tree_pairs[i]->update_map(0,part_taxa, false, true,central_node_part);
                
            } else {
                // for all partitions, which do not have the taxon of interest, remove the branch, which will be devided by the insertion of a new taxon from the backward map.
                // in the next step it will be substituted by the three branches: one with a new taxon and two ends of the devided branch
                for(j=0; j<nei_part_1->link_neighbors.size();j++){
                    if((nei_part_1->link_neighbors[j]->node->id == node_1_branch->id && nei_part_2->link_neighbors[j]->node->id == node_2_branch->id) or (nei_part_1->link_neighbors[j]->node->id == node_2_branch->id && nei_part_2->link_neighbors[j]->node->id == node_1_branch->id)){
                        
                        /*int end_v = nei_part_1->link_neighbors.size();
                        cout<<"======================================"<<endl;
                        cout<<"THE SUSPECT is at "<<j<<endl;
                        cout<<"BEFORE erase:"<<endl;
                        cout<<"nei_part_1->link_neighbors:"<<endl;
                        for(int h=0; h<nei_part_1->link_neighbors.size(); h++){
                            cout<<h<<":"<<nei_part_1->link_neighbors[h]<<endl;
                        }
                        cout<<"nei_part_2->link_neighbors:"<<endl;
                        for(int h=0; h<nei_part_2->link_neighbors.size(); h++){
                            cout<<h<<":"<<nei_part_2->link_neighbors[h]<<endl;
                        }
                        
                        nei_part_1->link_neighbors[j]=nullptr;
                        nei_part_2->link_neighbors[j]=nullptr;*/
                        
                        nei_part_1->link_neighbors.erase(nei_part_1->link_neighbors.begin()+j);
                        nei_part_2->link_neighbors.erase(nei_part_2->link_neighbors.begin()+j);
                        
                        /*cout<<"AFTER erase:"<<endl;
                        cout<<"nei_part_1->link_neighbors:"<<endl;
                        for(int h=0; h<nei_part_1->link_neighbors.size(); h++){
                            cout<<h<<":"<<nei_part_1->link_neighbors[h]<<endl;
                        }
                        cout<<"nei_part_2->link_neighbors:"<<endl;
                        for(int h=0; h<nei_part_2->link_neighbors.size(); h++){
                            cout<<h<<":"<<nei_part_2->link_neighbors[h]<<endl;
                        }
                        
                        
                        cout<<"THE SUSPECT:"<<endl;
                        cout<<j<<":"<<nei_part_1->link_neighbors[j]<<endl;
                        cout<<j<<":"<<nei_part_2->link_neighbors[j]<<endl;
                        cout<<"END of vector:"<<endl;
                        cout<<end_v-1<<":"<<nei_part_1->link_neighbors[end_v-1]<<endl;
                        cout<<end_v-1<<":"<<nei_part_2->link_neighbors[end_v-1]<<endl;
                        cout<<"Size of vectors = "<<nei_part_1->link_neighbors.size()<<" = "<<nei_part_2->link_neighbors.size()<<endl;
                        
                        //nei_part_1->link_neighbors.push_back(nullptr);
                        //nei_part_1->link_neighbors.push_back(nullptr);
                        
                        //cout<<"END of vector:"<<endl;
                        //cout<<end_v-1<<":"<<nei_part_1->link_neighbors[end_v-1]<<endl;
                        //cout<<end_v-1<<":"<<nei_part_2->link_neighbors[end_v-1]<<endl;
                        
                        //nei_part_1->link_neighbors.erase(nei_part_1->link_neighbors.begin()+nei_part_1->link_neighbors.size()-1);
                        //nei_part_2->link_neighbors.erase(nei_part_2->link_neighbors.begin()+nei_part_2->link_neighbors.size()-1);
                        
                         */
                         
                        break;
                    }
                }
            }
        } else if(part_tree_pairs[i]->findLeafName(node_name)){
            
            part_tree_pairs[i]->matrix->pr_ab_matrix[part_tree_pairs[i]->matrix->findTaxonID(node_name)][0]=1;
            
            if(induced_trees[i]->leafNum == 2){
                assert(induced_trees[i]->root->isLeaf() && "ERROR: in extendNewTaxon: root is not a leaf... something is wrong");
                induced_part_tree_branch_1[i]=induced_trees[i]->root;
                induced_part_tree_branch_2[i]=induced_trees[i]->root->neighbors[0]->node;
                induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) induced_part_tree_branch_1[i], (TerraceNode*) induced_part_tree_branch_2[i]);
                
                // Partitions with less than 3 taxa were not linked before. When you insert 3rd taxon, you should map them
                NodeVector part_taxa;
                part_taxa.clear();
                part_tree_pairs[i]->matrix->getPartTaxa(0, part_tree_pairs[i], part_tree_pairs[i]->induced_trees[0], part_taxa);
                part_tree_pairs[i]->linkTree(0, part_taxa, false, true);
                
            } else if(induced_trees[i]->leafNum == 1){
                TerraceNode* node_new = (TerraceNode*)induced_trees[i]->newNode(1, node_name.c_str());
                induced_trees[i]->root->addNeighbor(node_new, 0.0,0);
                node_new->addNeighbor(induced_trees[i]->root, 0.0,0);
                induced_trees[i]->leafNum += 1;
                induced_trees[i]->nodeNum += 1;
                induced_trees[i]->initializeTree();
                
            } else {
                induced_trees[i]->root = (TerraceNode*)induced_trees[i]->newNode(0, node_name.c_str());
                induced_trees[i]->leafNum = 1;
                induced_trees[i]->nodeNum = 1;
            }
        }
    }
    
    // DONE: I think, for corresponding branch, where a taxon has to be inserted, you need to set link_neighbors for all partitions to nullptr
    ((TerraceNeighbor*)node_1_branch->findNeighbor(node_2_branch))->link_neighbors.clear();
    ((TerraceNeighbor*)node_2_branch->findNeighbor(node_1_branch))->link_neighbors.clear();
    
    //cleanAllLinkNeighboursAndTaxa();
    
    insertNewTaxon(node_name,node_1_branch,node_2_branch);
    taxa_num += 1;
    matrix->extend_by_new_taxa(node_name, pr_ab_info);
    
    //assert(findLeafName(node_name) && "ERROR: Newly inserted leaf is not found! Insertion failed.");
    
    TerraceNode* leaf_node = (TerraceNode*) findLeafName(node_name);
    //cout<<"LEAF NAME_1:"<<leaf_node->name<<endl;
    TerraceNode* center_node = (TerraceNode*) leaf_node->neighbors[0]->node;
    
    //assert(center_node->findNeighbor(leaf_node));
    //cout<<"LEAF NAME_2:"<<leaf_node->name<<endl;
    //assert(leaf_node->name==node_name && "ERROR: Leaf name is incorrect.");

    //if(!findLeafName(node_name)){
    //    cout<<"ERROR: did not find the leaf node!"<<endl;
    //}
	if(!center_node->isNeighbor(leaf_node) or !leaf_node->isNeighbor(center_node) or !leaf_node->isLeaf()){
        printTaxa(cout);
		cout<<"Right after assignment:"<<endl;
		cout<<"LEAF:"<<leaf_node<<"|"<<leaf_node->id<<"|"<<leaf_node->name<<endl;
		cout<<"CENTER:"<<center_node<<"|"<<center_node->id<<"|"<<center_node->name<<endl;
    }
    
    TerraceNeighbor *center_node_nei;
    TerraceNeighbor *nei_aux;  
    // INFO: since you are introducing new branches, make sure the link_neighbor vector is initialised for them
    FOR_NEIGHBOR_IT(center_node, NULL, it){
        center_node_nei=(TerraceNeighbor*)(*it)->node->findNeighbor(center_node);
        center_node_nei->link_neighbors.resize(part_num,nullptr);
        ((TerraceNeighbor*)(*it))->link_neighbors.resize(part_num,nullptr);
    }

    //cout<<"Reached update part"<<endl;
    for(i=0; i<part_num; i++){
        if(induced_trees[i]->leafNum>2){
            if(induced_trees[i]->findLeafName(node_name)){
                //cout<<"Partition:"<<i<<"- larger than 2 - leaf occurs"<<endl;
                // if a taxon was inserted to the induced partition tree, update
                part_taxa.clear();
                matrix->getPartTaxa(i, this, induced_trees[i], part_taxa);
                update_map(i,part_taxa, true, false, center_node);

            }else{
                //cout<<"Partition:"<<i<<"- larger than 2 - leaf does not occur"<<endl;
                
                //cout<<"if a taxon does not occur on the induced partition tree"<<endl;
                nei_part_1 = (TerraceNeighbor*)induced_part_tree_branch_1[i]->findNeighbor(induced_part_tree_branch_2[i]);
                nei_part_2 = (TerraceNeighbor*)induced_part_tree_branch_2[i]->findNeighbor(induced_part_tree_branch_1[i]);
                
                /*//cout<<"First branch: one part of the devided branch"<<endl;
                nei_aux = (TerraceNeighbor*)node_1_branch->findNeighbor(center_node);
                nei_aux->link_neighbors[i] = nei_part_1;
                nei_part_1->link_neighbors.push_back(nei_aux);
                
                nei_aux = (TerraceNeighbor*)center_node->findNeighbor(node_1_branch);
                nei_aux->link_neighbors[i] = nei_part_2;
                nei_part_2->link_neighbors.push_back(nei_aux);
                
                //cout<<"Second branch: another part of the devided branch"<<endl;
                nei_aux = (TerraceNeighbor*)center_node->findNeighbor(node_2_branch);
                nei_aux->link_neighbors[i] = nei_part_1;
                nei_part_1->link_neighbors.push_back(nei_aux);
                
                nei_aux = (TerraceNeighbor*)node_2_branch->findNeighbor(center_node);
                nei_aux->link_neighbors[i] = nei_part_2;
                nei_part_2->link_neighbors.push_back(nei_aux);
                
                //cout<<"Third branch: incident to a newly inserted taxon: center->leaf"<<endl;
                nei_aux = (TerraceNeighbor*)center_node->findNeighbor(leaf_node);
                nei_aux->link_neighbors[i] = nei_part_1;
                nei_part_1->link_neighbors.push_back(nei_aux);
		
                //cout<<"Third branch: incident to a newly inserted taxon: leaf->center"<<endl;
                nei_aux = (TerraceNeighbor*)leaf_node->findNeighbor(center_node);
                nei_aux->link_neighbors[i] = nei_part_2;
                nei_part_2->link_neighbors.push_back(nei_aux);*/
                
            
                assert(center_node->neighbors.size()==3 && "ERROR: The central node does not have 3 neighbours! Case leafNum>2 & findLeafName(node_name) = false.");
                FOR_NEIGHBOR_IT(center_node, NULL, it){
                    
                    center_node_nei=(TerraceNeighbor*)(*it)->node->findNeighbor(center_node);

                    // backward map
                    nei_part_1->link_neighbors.push_back(center_node_nei);
                    nei_part_2->link_neighbors.push_back((TerraceNeighbor*)(*it));
                    
                    // forward map
                    center_node_nei->link_neighbors[i] = nei_part_1;
                    ((TerraceNeighbor*)(*it))->link_neighbors[i]=nei_part_2;
                }
            }
        }
        /*else{
            cout<<"Partition:"<<i<<"- less than 2"<<endl;
        }*/
    }
    
    
    //cout<<"INTERMEDIATE_INFO_TAXA_"<<leafNum<<"_INSERTED_"<<node_name<<"_TREE_"<<endl;
    //printTree(cout, WT_BR_SCALE | WT_NEWLINE);
    
    intermediated_trees_num +=1;
    if(intermediated_trees_num % 1000 == 0){
        cout<<"... trees generated - "<<intermediated_trees_num<<"; intermediated - "<<intermediated_trees_num-terrace_trees_num<<"; terrace - "<<terrace_trees_num<<"; dead paths - "<<dead_ends_num<<endl;
    }
    
    if(intermediated_trees_num - terrace_trees_num > intermediate_max_trees){
        write_warning_stop(1);
    }
    
    double time= getCPUTime()-Params::getInstance().startCPUTime;
    if(time > seconds_max){
        write_warning_stop(3);
    }
    
    
    // For DEBUGGING:
    /*bool clean_induced_part_maps = true;
    for(i=0; i<part_num; i++){
        part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
    }
    cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
    relinkALL(part_tree_pairs);*/
}

void Terrace::extendNewTaxon_naive(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> part_tree_pairs,vector<int> pr_ab_info){
    
    TerraceNeighbor *nei_1, *nei_2;
    nei_1 = (TerraceNeighbor*) node_1_branch->findNeighbor(node_2_branch);
    nei_2 = (TerraceNeighbor*) node_2_branch->findNeighbor(node_1_branch);
    
    Neighbor *nei_part_1, *nei_part_2;
    
    NodeVector part_taxa;
    
    int i;
    
    for(i=0; i<part_num; i++){
        if(part_tree_pairs[i]->findLeafName(node_name)){
            if(induced_trees[i]->leafNum>2){
                nei_part_1 = nei_1->link_neighbors[i];
                nei_part_2 = nei_2->link_neighbors[i];
                
                induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) nei_part_1->node, (TerraceNode*)nei_part_2->node);
                
            }else if(induced_trees[i]->leafNum == 2){
                assert(induced_trees[i]->root->isLeaf() && "ERROR: in extendNewTaxon: root is not a leaf... something is wrong");
                induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) induced_trees[i]->root, (TerraceNode*) induced_trees[i]->root->neighbors[0]->node);
            }else if(induced_trees[i]->leafNum == 1){
                TerraceNode* node_new = (TerraceNode*)induced_trees[i]->newNode(1, node_name.c_str());
                induced_trees[i]->root->addNeighbor(node_new, 0.0,0);
                node_new->addNeighbor(induced_trees[i]->root, 0.0,0);
                induced_trees[i]->leafNum += 1;
                induced_trees[i]->nodeNum += 1;
                induced_trees[i]->initializeTree();
            }else{
                assert(induced_trees[i]->leafNum == 0 && "ERROR in extendNewTaxon_naive: a tree should have 0 zero leaves...");
                induced_trees[i]->root = (TerraceNode*)induced_trees[i]->newNode(0, node_name.c_str());
                induced_trees[i]->leafNum += 1;
                induced_trees[i]->nodeNum += 1;
            }
            
            int id = part_tree_pairs[i]->matrix->findTaxonID(node_name);
            part_tree_pairs[i]->matrix->pr_ab_matrix[id][0] = 1;
            
        }
        part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(true);
    }
    
    cleanAllLinkNeighboursAndTaxa(true);
    insertNewTaxon(node_name, node_1_branch, node_2_branch);
    matrix->extend_by_new_taxa(node_name, pr_ab_info);
    taxa_num += 1;
    
    relinkALL(part_tree_pairs);
    
    intermediated_trees_num +=1;
    
}

void Terrace::generateTerraceTrees(Terrace *terrace, vector<Terrace*> part_tree_pairs, vector<string> *list_taxa_to_insert, int taxon_to_insert){
    string taxon_name;
    taxon_name = list_taxa_to_insert->at(taxon_to_insert);
    //cout<<endl<<"*******************************************************"<<endl;
    //cout<<"GENERATE_TERRACE_TREES | TAXON "<<taxon_name<<endl;
    //cout<<"*******************************************************"<<endl;
    
    NodeVector node1_vec_branch, node2_vec_branch;
    int j, id;
    getAllowedBranches(taxon_name, part_tree_pairs, &node1_vec_branch, &node2_vec_branch);
    
    if(!node1_vec_branch.empty()){
        //cout<<"NUM_OF_ALLOWED_BRANCHES_"<<taxon_name<<"_"<<node1_vec_branch.size()<<endl;
        //cout<<"ALL ALLOWED BRANCHES:"<<endl;
        //for(j=0; j<node1_vec_branch.size(); j++){
        //    cout<<j<<":"<<node1_vec_branch[j]->id<<"-"<<node2_vec_branch[j]->id<<endl;
        //}
        
        for(j=0; j<node1_vec_branch.size(); j++){
            //cout<<"-----------------------------------"<<endl<<"INSERTing taxon "<<taxon_name<<" on branch "<<j+1<<" out of "<<node1_vec_branch.size()<<": "<<node1_vec_branch[j]->id<<"-"<<node2_vec_branch[j]->id<<endl<<"-----------------------------------"<<endl;
            
            id = terrace->matrix->findTaxonID(taxon_name);
            assert(id!=-1);
            extendNewTaxon(taxon_name,(TerraceNode*)node1_vec_branch[j],(TerraceNode*)node2_vec_branch[j],part_tree_pairs,terrace->matrix->pr_ab_matrix[id]);
            //extendNewTaxon_naive(taxon_name,(TerraceNode*)node1_vec_branch[j],(TerraceNode*)node2_vec_branch[j],part_tree_pairs,terrace->matrix->pr_ab_matrix[id]);
            
            if(taxon_to_insert != list_taxa_to_insert->size()-1){
                
                generateTerraceTrees(terrace, part_tree_pairs, list_taxa_to_insert, taxon_to_insert+1);
                
                // INFO: IF NEXT TAXON DOES NOT HAVE ALLOWED BRANCHES CURRENT TAXON IS DELETED AND NEXT BRANCH IS EXPLORED.
                //remove_one_taxon_naive(taxon_name,part_tree_pairs);
                remove_one_taxon(taxon_name,part_tree_pairs);
                
            } else {
                if(terrace_out){
                    if(trees_out_lim==0 or terrace_trees.size()<trees_out_lim){
                    terrace_trees.push_back(getTreeTopologyString(this));
                    //ofstream out;
                    //out.exceptions(ios::failbit | ios::badbit);
                    //out.open(out_file,std::ios_base::app);
                    //printTree(out, WT_BR_SCALE | WT_NEWLINE);
                    //out.close();
                    }
                }
                terrace_trees_num+=1;
                //if(terrace_trees_num % 1000 == 0){
                //    cout<<"... generated tree "<<terrace_trees_num<<endl;
                //}
                //printTree(cout, WT_BR_SCALE | WT_NEWLINE);
                if(terrace_trees_num == terrace_max_trees){
                    write_warning_stop(2);
                }
                //remove_one_taxon_naive(taxon_name,part_tree_pairs);
                remove_one_taxon(taxon_name,part_tree_pairs);
            }
        }
    } else {
        //cout<<"NUM_OF_ALLOWED_BRANCHES_"<<taxon_name<<"_0_dead_end"<<endl;
        //cout<<"For a given taxon "<<taxon_name<<" there are no allowed branches.. Dead end.."<<endl;
        dead_ends_num +=1;
    }

}

void Terrace::remove_one_taxon(string taxon_name, vector<Terrace*> part_tree_pairs){
    
    //cout<<"-----------------------------------"<<endl<<"REMOVING TAXON: "<<taxon_name<<endl<<"-----------------------------------"<<endl;
    
    int i, j, h;
    NodeVector induced_part_tree_branch_1, induced_part_tree_branch_2;
    induced_part_tree_branch_1.resize(part_num,nullptr);
    induced_part_tree_branch_2.resize(part_num,nullptr);
    
    Node *node_1, *node_2; // branch nodes of the agile (parent) tree, that will be the ends of the joint branch after taxon removal
    TerraceNeighbor *neiC1, *nei1C;
    
    Node *node_leaf_main, *node_central_main;
    Node *node_leaf, *node_central;
    NodeVector branch_nodes;
    TerraceNeighbor *nei1, *nei2, *nei_aux;
    
    node_leaf_main = findLeafName(taxon_name);
    node_central_main = node_leaf_main->neighbors[0]->node;
    FOR_NEIGHBOR_DECLARE(node_central_main, node_leaf_main, it){
        branch_nodes.push_back((*it)->node);
    }
    node_1 = branch_nodes[0];
    node_2 = branch_nodes[1];
    branch_nodes.clear();
    
    // one of the two branches that will be joint on the agile tree
    neiC1 = (TerraceNeighbor*) node_central_main->findNeighbor(node_1);
    nei1C = (TerraceNeighbor*) node_1->findNeighbor(node_central_main);
    
    // Below NeiVectors will store neibours from agile tree and top induced trees that need to be updated
    // Corresponding neibours will have their forward map equal to nei_low_12 on corresponding low induced subtree, i.e. I'm going to "fold" all nei's in the same direction
    // nei_12 ->link = nei_low_12, while nei_21 ->link = nei_low_21
    vector<NeighborVec> update_nei; // neis on agile tree, for each partition its own set of neis to be updated
    update_nei.resize(part_num);
    
    vector<NeighborVec> update_nei_top; // neis on top induced partition tree, for each pair its own set of neis to be updated
    update_nei_top.resize(part_num);
    
    NeighborVec update_nei_aux, update_nei_top_aux;
    bool appear;
    
    for(i=0; i<part_num; i++){
        //cout<<"-------------------------------"<<endl;
        //cout<<endl<<"Partition "<<i<<":"<<endl;
        
        update_nei_aux.clear();
        update_nei_top_aux.clear();
        
        //part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(true);
        
        /*
         *  If the taxon to be removed is present in partition, clear corresponding maps (forward and backward maps),
         *  which involve three branches (6 neighbors): incident to the taxon and two parts of the branch to be joint after taxon removal.
         *  For all partitions collect branch nodes of common partition subtrees, whose maps (forward or backward) have to
         *  be modified (clear and update).
         */
        
        if(part_tree_pairs[i]->findLeafName(taxon_name)){
            //cout<<"- with Leaf"<<endl;
            
            int id = part_tree_pairs[i]->matrix->findTaxonID(taxon_name);
            part_tree_pairs[i]->matrix->pr_ab_matrix[id][0] = 0;
            
            if(induced_trees[i]->leafNum>2){
                
                //cout<<"- with more than 2 leaves:"<<induced_trees[i]->leafNum<<endl;
                node_leaf = induced_trees[i]->findLeafName(taxon_name);
                node_central = node_leaf->neighbors[0]->node;
                
                
                branch_nodes.clear();
                FOR_NEIGHBOR(node_central, node_leaf, it){
                    branch_nodes.push_back((*it)->node);
                }
                
                induced_part_tree_branch_1[i] = branch_nodes[0];
                induced_part_tree_branch_2[i] = branch_nodes[1];
                
                if(neiC1->link_neighbors[i]->node == branch_nodes[1] or nei1C->link_neighbors[i]->node == branch_nodes[1]){
                    induced_part_tree_branch_1[i] = branch_nodes[1];
                    induced_part_tree_branch_2[i] = branch_nodes[0];
                }
                
                /*if(neiC1->link_neighbors[i]->id == neiC2->link_neighbors[i]->id){
                    induced_part_tree_branch_1[i] = neiC1->link_neighbors[i]->node;
                    induced_part_tree_branch_2[i] = ((TerraceNeighbor*)node_1->findNeighbor(node_central_main))->link_neighbors[i]->node;
                    
                    if(induced_part_tree_branch_1[i] == node_leaf or induced_part_tree_branch_2[i] == node_leaf){
                        branch_nodes.clear();
                        FOR_NEIGHBOR(node_central, node_leaf, it){
                            branch_nodes.push_back((*it)->node);
                        }
                        induced_part_tree_branch_1[i] = branch_nodes[0];
                        induced_part_tree_branch_2[i] = branch_nodes[1];
                    }
                }else{
                    induced_part_tree_branch_1[i] = neiC1->link_neighbors[i]->node;
                    induced_part_tree_branch_2[i] = neiC2->link_neighbors[i]->node;
                }*/
                
                //cout<<"induced_part_tree_branch_1[i] = "<<induced_part_tree_branch_1[i]->id<<endl;
                //cout<<"induced_part_tree_branch_2[i] = "<<induced_part_tree_branch_2[i]->id<<endl;
                //cout<<"node_leaf="<<node_leaf->id<<endl;
                //cout<<"node_central="<<node_central->id<<endl;
                
                FOR_NEIGHBOR(node_central, nullptr, it){
                    nei1=(TerraceNeighbor*)(*it);
                    nei2=(TerraceNeighbor*)(*it)->node->findNeighbor(node_central);
                    
                    if(!((TerraceNeighbor*)(*it))->link_neighbors.empty()){
                        //cout<<endl<<"DEBUGGING: clearing forward nei's maps"<<endl;
                        // clearing forward maps from agile (parent) tree to common partition subtree
                        for(j=0; j<nei1->link_neighbors.size(); j++){
                            //cout<<"Neibour "<<j<<":"<<endl;
                            //cout<<"- nei1->link_neighbors[j]="<<nei1->link_neighbors[j]->node->id<<endl;
                            //cout<<"- nei2->link_neighbors[j]="<<nei2->link_neighbors[j]->node->id<<endl;
                            ((TerraceNeighbor*)nei1->link_neighbors[j])->link_neighbors[i] = nullptr;
                            ((TerraceNeighbor*)nei2->link_neighbors[j])->link_neighbors[i] = nullptr;
                            
                            
                            // INFO: you do not want any of the 6 neibours that will be modified to appear in the list of neis to be updated, because these neis won't exist after removal and you'll get SEG_FAULT
                            appear=false;
                            FOR_NEIGHBOR_DECLARE(node_central_main, nullptr, it1){
                                if(nei1->link_neighbors[j] == (*it1)){
                                    appear = true;
                                    break;
                                }
                                if(nei1->link_neighbors[j] == (*it1)->node->findNeighbor(node_central_main)){
                                    appear = true;
                                    break;
                                }
                            }
                            if(!appear){
                                update_nei_aux.push_back(nei1->link_neighbors[j]);
                                //cout<<"+ added nei1->"<<nei1->link_neighbors[j]->node->id<<" to update list"<<endl;
                            }
                            
                            appear=false;
                            FOR_NEIGHBOR(node_central_main, nullptr, it1){
                                if(nei2->link_neighbors[j] == (*it1)){
                                    appear = true;
                                    break;
                                }
                                if(nei2->link_neighbors[j] == (*it1)->node->findNeighbor(node_central_main)){
                                    appear = true;
                                    break;
                                }
                            }
                            if(!appear){
                                update_nei_aux.push_back(nei2->link_neighbors[j]);
                                //cout<<"+ added nei2->"<<nei2->link_neighbors[j]->node->id<<" to update list"<<endl;
                            }
                        }
                        
                        // clearing backward maps from common partition subtree to agile (parent) tree
                        nei1->link_neighbors.clear();
                        nei2->link_neighbors.clear();
                    }
                
                    if(!((TerraceNeighbor*)(*it))->link_neighbors_lowtop_back.empty()){
                        
                        //cout<<"clearing forward maps from induced partition tree (parent) to common partition subtree"<<endl;
                        for(j=0; j<((TerraceNeighbor*)(*it))->link_neighbors_lowtop_back.size(); j++){
                            ((TerraceNeighbor*)nei1->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                            ((TerraceNeighbor*)nei2->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                            
                            // INFO: note, that for top partition trees we do not need to avoid adding 6 neis to the list, because they will still exist after taxon deletion. However, when you'll update, set all neis to nei_low_12 and nei_low_21 correctly to nei_low_21 and not to nei_low_12
                            update_nei_top_aux.push_back(nei1->link_neighbors_lowtop_back[j]);
                            update_nei_top_aux.push_back(nei2->link_neighbors_lowtop_back[j]);
                            //cout<<"update_nei_top_aux.size()="<<update_nei_top_aux.size()<<endl;
                        }
                        
                        // clearing backward maps from common partition subtree to induced partition tree (parent)
                        nei1->link_neighbors_lowtop_back.clear();
                        nei2->link_neighbors_lowtop_back.clear();
                    }
                }
                
                
                // Get your "to be updated list of neighbours" in place for agile tree
                update_nei[i] = update_nei_aux;
                
                // Get your "to be updated list of neighbours" in place for top induced partition tree
                update_nei_top[i] = update_nei_top_aux;
                
            } //else {
              //  cout<<"- with less than 2 leaves:"<<induced_trees[i]->leafNum<<endl;
            //}
        } else {
            //cout<<"- without Leaf"<<endl;
            if(induced_trees[i]->leafNum>2){
                //cout<<"- with more than 2 leaves:"<<induced_trees[i]->leafNum<<endl;
                
                // CHECK: again, this might by the same node... in case things are messed up
                induced_part_tree_branch_1[i] = neiC1->link_neighbors[i]->node;
                //induced_part_tree_branch_2[i] = neiC2->link_neighbors[i]->node;
                induced_part_tree_branch_2[i] = nei1C->link_neighbors[i]->node;
                
                //cout<<"induced_part_tree_branch_1[i] = "<<induced_part_tree_branch_1[i]->id<<endl;
                //cout<<"induced_part_tree_branch_2[i] = "<<induced_part_tree_branch_2[i]->id<<endl;
                
                /*
                 * INFO: For partitions, which do not have corresponding taxon, remove 6 possible backward neighbors corresponding to three branches, which will be removed/modified
                 */
                nei1 = (TerraceNeighbor*)neiC1->link_neighbors[i];
                nei2 = (TerraceNeighbor*)nei1C->link_neighbors[i];
                
                for(j=nei1->link_neighbors.size()-1; j!=-1;j--){
                    for(h=0; h<node_central_main->neighbors.size();h++){
                        nei_aux = (TerraceNeighbor*)(node_central_main->neighbors[h]->node->findNeighbor(node_central_main));
                        if(nei1->link_neighbors[j] == node_central_main->neighbors[h]){
                            nei1->link_neighbors.erase(nei1->link_neighbors.begin()+j);
                            assert(nei2->link_neighbors[j] == nei_aux);
                            nei2->link_neighbors.erase(nei2->link_neighbors.begin()+j);
                        } else if(nei2->link_neighbors[j] == node_central_main->neighbors[h]){
                            nei2->link_neighbors.erase(nei2->link_neighbors.begin()+j);
                            assert(nei1->link_neighbors[j] == nei_aux);
                            nei1->link_neighbors.erase(nei1->link_neighbors.begin()+j);
                        }
                    }
                }
            }
            //else{
            //    cout<<"- with less than 2 leaves:"<<induced_trees[i]->leafNum<<endl;
            //}
        }
    }
    
    //cleanAllLinkNeighboursAndTaxa(true);

    //cout<<"Removing the taxon from the low-level induced partition tree...."<<endl;
    for(i=0; i<part_num; i++){
        if(induced_trees[i]->leafNum>1){
            if(induced_trees[i]->findLeafName(taxon_name)){
                induced_trees[i]->remove_taxon(taxon_name);
            }
        } else if(induced_trees[i]->root){
            if(induced_trees[i]->root->name == taxon_name){
                induced_trees[i]->remove_taxon(taxon_name);
            }
        }
    }

    //cout<<"Removing the taxon from agile tree...."<<endl;
    remove_taxon(taxon_name);
    
    ((TerraceNeighbor*)node_1->findNeighbor(node_2))->link_neighbors.resize(part_num,nullptr);
    ((TerraceNeighbor*)node_2->findNeighbor(node_1))->link_neighbors.resize(part_num,nullptr);
    
    matrix->remove_taxon(taxon_name);
    taxa_num -= 1;
    
    //cout<<"Updating maps for partition trees...."<<endl;
    // RE-MAP with new UPDATE using vectors of "to be modified" and not the update_map function
    for(i=0; i<part_num; i++){
        //cout<<"Updating maps for partition "<<i<<":"<<endl;
        if(induced_trees[i]->leafNum>2){
            nei1 = (TerraceNeighbor*) induced_part_tree_branch_1[i]->findNeighbor(induced_part_tree_branch_2[i]);
            nei2 = (TerraceNeighbor*) induced_part_tree_branch_2[i]->findNeighbor(induced_part_tree_branch_1[i]);
            
            // AGILE TREE
            if(!update_nei[i].empty()){
                //cout<<"There are "<<update_nei[i].size()/2<<" branches ("<<update_nei[i].size()<<" neighbors) to be updated"<<endl;
                for(j=0;j<update_nei[i].size()/2; j++){
                    ((TerraceNeighbor*)update_nei[i][2*j])->link_neighbors[i] = nei1;
                    nei1->link_neighbors.push_back(update_nei[i][2*j]);
                    ((TerraceNeighbor*)update_nei[i][2*j+1])->link_neighbors[i] = nei2;
                    nei2->link_neighbors.push_back(update_nei[i][2*j+1]);
                }
            }
            ((TerraceNeighbor*)node_1->findNeighbor(node_2))->link_neighbors[i] = nei1;
            nei1->link_neighbors.push_back(node_1->findNeighbor(node_2));
            ((TerraceNeighbor*)node_2->findNeighbor(node_1))->link_neighbors[i] = nei2;
            nei2->link_neighbors.push_back(node_2->findNeighbor(node_1));
            
            // TOP induced partition trees
            if(!update_nei_top[i].empty()){
                //cout<<"TOP PART TREE: There are "<<update_nei_top[i].size()/2<<" branches ("<<update_nei_top[i].size()<<" neighbors) to be updated"<<endl;
                for(j=0;j<update_nei_top[i].size()/2; j++){
                    ((TerraceNeighbor*)update_nei_top[i][2*j])->link_neighbors[0] = nei1;
                    nei1->link_neighbors_lowtop_back.push_back(update_nei_top[i][2*j]);
                    ((TerraceNeighbor*)update_nei_top[i][2*j+1])->link_neighbors[0] = nei2;
                    nei2->link_neighbors_lowtop_back.push_back(update_nei_top[i][2*j+1]);
                }
            }
        }
    }
}

void Terrace::remove_one_taxon_naive(string taxon_name, vector<Terrace*> part_tree_pairs){
    //cout<<"-----------------------------------"<<endl<<"REMOVING TAXON: "<<taxon_name<<endl<<"-----------------------------------"<<endl;
    
    bool clean_induced_part_maps = true;
    
    for(int i=0; i<part_num; i++){
        part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
        if(part_tree_pairs[i]->findLeafName(taxon_name)){
            int id = part_tree_pairs[i]->matrix->findTaxonID(taxon_name);
            part_tree_pairs[i]->matrix->pr_ab_matrix[id][0] = 0;
        }
    }
    
    cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
    
    for(int i=0; i<part_num; i++){
        if(induced_trees[i]->findLeafName(taxon_name)){
            induced_trees[i]->remove_taxon(taxon_name);
        }
    }
    
    remove_taxon(taxon_name);
    matrix->remove_taxon(taxon_name);
    taxa_num -= 1;
    
    relinkALL(part_tree_pairs);
}

void Terrace::relinkALL(vector<Terrace*> part_tree_pairs){
    
    // re-link
    linkTrees(true, false);
    
    for(int i=0; i<part_tree_pairs.size(); i++){
        part_tree_pairs[i]->linkTrees(false, true);
    }
}

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
    
    if(node->empty_br_node_nei.size()>0){
        //node->empty_branches.clear();
        node->empty_br_node_nei.clear();
        node->empty_br_dad_nei.clear();
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        clearEmptyBranchAndTaxaINFO((TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
}

void Terrace::cleanAllLinkNeighboursAndTaxa(bool clean_induced_part_maps){
    cleanAllLinkINFO(clean_induced_part_maps);
}

void Terrace::print_ALL_DATA(vector<Terrace*> part_tree_pairs){
    int i;
    
    cout<<endl<<"================ BEGIN: PRINTING all INFO ================"<<endl;
    cout<<"Initial PARENT tree"<<endl;
    print_terrace_tree();
    
    matrix->print_pr_ab_matrix();
    
    cout<<"================ INDUCED PARTITION PAIRS =========="<<endl;
    for(i=0; i<part_num; i++){
        cout<<"------------------> TOP PART TREE "<<i<<":"<<endl;
        part_tree_pairs[i]->print_terrace_tree();
        part_tree_pairs[i]->matrix->print_pr_ab_matrix();
        
        cout<<endl<<"------------------> LOW part tree "<<i<<":"<<endl;
        induced_trees[i]->print_terrace_tree();
        cout<<endl;
        cout<<endl<<"========================================================="<<endl;
    }
    
    printMapInfo();
    printBackMapInfo();
    cout<<"================ END: PRINTING all INFO ================="<<endl<<endl;

}

void Terrace::renameTaxa(){
    
    getTaxaName(taxa_names_orgn);
    
    NodeVector taxa;
    getTaxa(taxa);
    
    for(int i=0; i<taxa_num; i++){
        stringstream ss;
        ss << i;
        string name = "t" + ss.str();
        for(NodeVector::iterator it = taxa.begin(); it!=taxa.end(); it++){
            if((*it)->name == taxa_names_orgn[i]){
                (*it)->name = name;
                break;
            }
        }
        for(int j=0; j<taxa_num; j++){
            if(matrix->taxa_names[j] == taxa_names_orgn[i]){
                matrix->taxa_names[j] = name;
                break;
            }
        }
    }
}

bool Terrace::check_two_trees(MTree* query_tree){
    TerraceTree tree;
    tree.copyTree_byTaxonNames(query_tree, matrix->taxa_names);
    Terrace *query_terrace = new Terrace(tree, matrix);
    
    vector<string> trees, taxa_names;
    
    for(int part=0; part<part_num; part++){
        
        trees.clear();
        trees.push_back(getTreeTopologyString(induced_trees[part]));
        //trees.push_back(getTreeTopologyString(induced_trees_query[part]));
        
        taxa_names.clear();
        induced_trees[part]->getTaxaName(taxa_names);
        
        MTreeSet tree_set;
        tree_set.init(trees,taxa_names,induced_trees[part]->rooted);
        //tree_set[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
        
        MTreeSet tree_set_1;
        trees.clear();
        trees.push_back(getTreeTopologyString(query_terrace->induced_trees[part]));
        tree_set_1.init(trees,taxa_names,query_terrace->induced_trees[part]->rooted);
        //tree_set_1[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
        
        int *rfdist;
        int n=1;
        rfdist = new int [n*n];
        memset(rfdist, 0, n*n* sizeof(int));
        tree_set.computeRFDist(rfdist,&tree_set_1);
        //cout<<"Partition "<<part+1<<": RF-distance between induced trees is "<<rfdist[0]<<endl;
        
        if(rfdist[0]!=0){
            //cout<<"---------------------------------------------"<<endl;
            //cout<<"A query tree does not belong to the terrace:"<<endl;
            //tree->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            
            //cout<<"The first encountered pair of induced partition trees that do not match. Partition "<<part+1<<endl;
            //tree_set[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            //tree_set_1[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            
            delete [] rfdist;
            tree_set.clear();
            return false;
        }
        delete [] rfdist;
        tree_set.clear();
        
    }
    //cout<<"---------------------------------------------"<<endl;
    //cout<<"A query tree is ON the terrace"<<endl;
    
    return true;
}

void Terrace::write_terrace_trees_to_file(){
 
    cout<<"Wall-clock time used so far: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")"<<endl;
    cout<<"CPU time used on generation of terrace trees: "
    << getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")" << endl;
    cout<<"---------------------------------------------------------"<<endl;
    
    cout<<"Printing "<<terrace_trees_num<<" terrace trees to file "<<endl<<out_file<<"..."<<endl;
        
    if(trees_out_lim==0 or trees_out_lim > terrace_trees_num){
        trees_out_lim = terrace_trees_num;
    } else {
        cout<<"WARNING: The number of generated trees from the terrace ("<<terrace_trees_num<<") is larger than the output treshold ("<<trees_out_lim<<" trees). Only "<<trees_out_lim<<" trees will be written to the file."<<endl;
    }
    
    ofstream out;
    out.exceptions(ios::failbit | ios::badbit);
    out.open(out_file,std::ios_base::app);
        
    for(int i=0; i<trees_out_lim; i++){
        out<<terrace_trees[i]<<endl;
    }
    out.close();

}

void Terrace::write_summary_generation(){

    Params::getInstance().run_time = (getCPUTime() - Params::getInstance().startCPUTime);
    
    cout<<"---------------------------------------------------------"<<endl;
    cout<<"SUMMARY:"<<endl;
    cout<<"Number of trees on terrace: "<<terrace_trees_num<<endl;
    cout<<"Number of intermediated trees visited: "<<intermediated_trees_num - terrace_trees_num<<endl;
    cout<<"Number of dead ends encountered: "<<dead_ends_num<<endl;
    cout<<"---------------------------------------------------------"<<endl;
    
    if(terrace_out){
        write_terrace_trees_to_file();
        cout<<"---------------------------------------------------------"<<endl;
    }
    cout<<"Total wall-clock time used: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")"<<endl;
    cout<<"Total CPU time used: "
    << getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")" << endl;
    cout<<endl;
}


void Terrace::write_warning_stop(int type){
    
    cout<<endl<<"=========================================="<<endl;
    cout<<"WARNING: stopping condition is active!"<<endl;
    cout<<"The total number of trees on the terrace is NOT yet computed!"<<endl;
    cout<<"Check summary at the current step. You can change the stopping rule via corresponding option (see below)."<<endl;
    cout<<"=========================================="<<endl;
    
    switch (type) {
        case 1:
            cout<<"Type of stopping rule: number of visited intermediate trees"<<endl;
            cout<<"Current setting: stop if the number of visited intermediate trees reached "<<intermediate_max_trees<<endl;
            cout<<"To change the value use -t_stop_i <number_of_intermediate_trees_to_stop>"<<endl;
            cout<<"------------------------------------------"<<endl;
            cout<<"The number of considered intermediated trees already reached "<<intermediate_max_trees<<". Exiting generation process..."<<endl;
            break;
        case 2:
            cout<<"Type of stopping rule: terrace size"<<endl;
            cout<<"Current setting: stop if the number of trees on the terrace reached "<<terrace_max_trees<<endl;
            cout<<"To change the value use -t_stop_t <number_of_terrace_trees_to_stop>"<<endl;
            cout<<"------------------------------------------"<<endl;
            cout<<"Considered terrace already contains "<<terrace_max_trees<<" trees."<<endl<<"Exiting generation process..."<<endl;
            break;
            
        case 3:
            cout<<"Type of stopping rule: CPU time used"<<endl;
            cout<<"Current setting: stop if the CPU time used is larger than "<<seconds_max<<endl;
            cout<<"To change the value use -t_stop_h <number_of_hours_to_stop>"<<endl;
            cout<<"------------------------------------------"<<endl;
            cout<<"Total CPU time is already larger than "<<seconds_max<<" trees."<<endl<<"Exiting generation process..."<<endl;
            break;
            
        default:
            break;
    }
    
    cout<<"Printing summary at current step..."<<endl;
    write_summary_generation();
    time_t end_time;
    time(&end_time);
    cout << "Date and Time: " << ctime(&end_time);
    exit(0);
        
}
