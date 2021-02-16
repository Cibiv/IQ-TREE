//
//  terracetree.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "terracetree.hpp"
#include "terracenode.hpp"

TerraceTree::TerraceTree(): MTree(){};
TerraceTree::~TerraceTree(){
    if (root != NULL)
        freeNode();
    root = NULL;
};


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
    //bool check;
    
    string taxa_set = "";
    NodeVector taxa_nodes;
    
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(tree->leafNum,0);
    
    tree->getTaxa(taxa_nodes);
    
    for(j=0; j<taxa_nodes.size(); j++){
        //check = false;
        for(i=0; i<taxa_num; i++){
            if(taxa_nodes[j]->name == taxon_names[i]){
                check_int[taxa_nodes[j]->id]=1;
                sum+=1;
                //check = true;
                break;
            }
        }
        //cout<<"Taxon["<<j<<"] = "<<check_int[j]<<"| main_tree:"<<taxa_nodes[j]->name<<"("<<taxa_nodes[j]->id<<") "<<"-> subtree: "<<((check) ? taxon_names[i] : "")<<endl;
    }
    assert(sum == taxa_num && "Not all of the taxa appear in the complete tree!");
    taxa_set.clear();
    taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
    copyTree(tree,taxa_set);

    //printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
}

void TerraceTree::cleanAllLinkINFO(bool clean_induced_part_maps, TerraceNode *node, TerraceNode *dad){

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
            //cout<<"| IF link_neighbours_exist -> clear them: size "<<nei->link_neighbors.size()<<endl;
            
            // Clearing backward map from induced partition trees...
            if(clean_induced_part_maps){
                for(part=0; part<nei->link_neighbors.size(); part++){
                    // INFO: since for the trees with less than 3 taxa you do not do any maps, first check if there is a link_neighbor for the neighbor on the parent tree
                    if((TerraceNeighbor*)nei->link_neighbors[part]){
                        if(((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.size()>0){
                            ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.clear();
                            ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors.clear();
                        }
                        if(((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors_lowtop_back.size()>0){
                            ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors_lowtop_back.clear();
                            ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors_lowtop_back.clear();
                        }
                        
                        /*if(((TerraceNeighbor*)nei->link_neighbors[part])->taxa_to_insert.size()>0){
                         ((TerraceNeighbor*)nei->link_neighbors[part])->taxa_to_insert.clear();
                         ((TerraceNeighbor*)dad_nei->link_neighbors[part])->taxa_to_insert.clear();
                         }*/
                    }
                }
            }
            
            // Clearing forward map from parent tree to induced partition trees...
            nei->link_neighbors.clear();
            dad_nei->link_neighbors.clear();
            
        }
        
        //node->empty_taxa.clear();
        //node->empty_branches.clear();
        node->empty_br_dad_nei.clear();
        node->empty_br_node_nei.clear();
        
        //dad->empty_taxa.clear();
        //dad->empty_branches.clear();
        dad->empty_br_dad_nei.clear();
        dad->empty_br_node_nei.clear();
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it){
        cleanAllLinkINFO(clean_induced_part_maps, (TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
}

void TerraceTree::insertNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch){
    
    // IMPORTANT!!!!!! WARNING: CAREFULL about calling a presence_absence matrix entries by the leaf id!!!!!! CHECK the corresponding code, because your maps will be fucked up!
    // TODO: safe option is a function to re-number taxa, nodes -> I'm using the initializeTree();
    
    TerraceNode *node_1, *node_2;
    
    node_1 = (TerraceNode*) newNode(leafNum, node_name.c_str());
    
    leafNum += 1;
    nodeNum += 1;
    
    node_2 = (TerraceNode*)(newNode(nodeNum));
    nodeNum += 1;
    
    int br_id = branchNum;
    branchNum += 1;
    
    // WARNING: you need unique branch ids for intersection in allowed branches!!!!!
    // CHECK: if there are any issues with current handling of branch ids -> using initializeTree();
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
    
    initializeTree();
    
    //drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    
}

void TerraceTree::remove_taxon(string taxon_name){
    
    // WARNING: I think it does not remove taxon correctly and some neighbours appear to be present as Neighbors instead of TerraceNeighbours

    StrVector taxa;
    taxa.push_back(taxon_name);
  
    if(leafNum>=2){
        TerraceNode *node = (TerraceNode*) findLeafName(taxon_name);
        TerraceNeighbor* nei = (TerraceNeighbor*)node->neighbors[0];

        // nei is a central node (or a second leaf in a 2-taxon tree), the branch will be joined, therefore, the link_neis should be deleted (only the pointers, not objects).
        FOR_NEIGHBOR_DECLARE(nei->node,NULL, it){
            ((TerraceNeighbor*)(*it))->delete_ptr_members();
            ((TerraceNeighbor*)(*it)->node->findNeighbor(nei->node))->delete_ptr_members();
        }
    }
        
    if(leafNum>2){
        removeTaxa(taxa);
        initializeTree();
    } else {
        if(leafNum==2){
            //cout<<"two-taxon tree, remove one taxon"<<endl;

            if(root->name == taxon_name){
                TerraceNode * new_root = (TerraceNode*) root->neighbors[0]->node;
                delete root;
                
                // Free pointers for the remaining node:
                new_root->deleteNode();
                
                // Set new root
                root = new_root;
                leafNum = 1;
                nodeNum = 1;
                branchNum = 0;
                root->id = 0;
            }else{
                delete root->neighbors[0]->node;
                
                // Free pointers for the remaining node:
                root->deleteNode();
                
                leafNum = 1;
                nodeNum = 1;
                branchNum = 0;
                root->id = 0;
            }
        }else if(leafNum==1){
            delete root;
            root = nullptr;
            leafNum = 0;
            nodeNum = 0;
            branchNum = 0;
        }
    }

}

void TerraceTree::print_terrace_tree(bool draw){
    
    if(leafNum>2 && draw){
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    }else if(leafNum==2 or (!draw && leafNum>1)){
        printTree(cout, WT_BR_SCALE | WT_NEWLINE);
    }else if(leafNum==1){
        cout<<"("<<root->name<<");"<<endl;
    }else{
        cout<<"();"<<endl;
    }
    
}

string getTreeTopologyString(MTree* tree){
    stringstream tree_stream;
    tree->printTree(tree_stream, WT_BR_LEN_ROUNDING + WT_SORT_TAXA);
    return tree_stream.str();
}
