//
//  presenceabsencematrix.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "presenceabsencematrix.hpp"
#include "tree/node.h"
#include "tree/mtree.h"


void PresenceAbsenceMatrix::read_pr_ab_matrix(const char *infile){
    ifstream in;
    //cout<<endl<<"-----------------------------------------------------"<<endl;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        read_pr_ab_matrix(in);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
}

void PresenceAbsenceMatrix::read_pr_ab_matrix(istream &in){
    string str_rest,name;
    
    if(!(in>>taxa_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    if(!(in>>part_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    
    int i=0,j=0;
    for(i=0; i<taxa_num; i++){
        if(!(in>>name)) throw "Each line should start with a taxon name!";
        if(name == "0" or name == "1") throw "Each line should start with a taxon name! 0 and 1 are not allowed as taxon names.";
        taxa_names.push_back(name);
        vector<int> vec(part_num, -1);
        for(j=0; j<part_num; j++){
            if(!(in >> vec[j])) throw "Could not read a matrix entry! For each species make sure there are as many entries as the number of partitions specified in the first line of the file. Moreover, presence-absence matrix should only contain 0, 1!";
            if(vec[j] < 0) throw "Error: A negative entry! Presence-absence matrix should only contain 0, 1!";
            if(vec[j] > 1) throw "Error: The entry is greater than 1! Presence-absence matrix should only contain 0, 1!";
        }
        pr_ab_matrix.push_back(vec);
    }
    
    init();
};

void PresenceAbsenceMatrix::print_pr_ab_matrix(){
    
    cout<<"Presence-absence matrix:"<<endl;
    for(int i=0; i<taxa_num; i++){
        cout<<taxa_names[i]<<" ";
        for(int j=0; j<part_num; j++){
            cout<<pr_ab_matrix[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
};


int PresenceAbsenceMatrix::findTaxonID(string taxon_name){
    int id;
    for(id=0; id<taxa_names.size(); id++){
        if(taxa_names[id]==taxon_name){
            //cout<<"GET_TAXON_ID: taxa_name = "<<taxon_name<<" | taxa_id = "<<id<<endl;
            return id;
        }
    }
    return -1;
};

void PresenceAbsenceMatrix::init(){
    flag_reorderAccordingToTree = false;
}

void PresenceAbsenceMatrix::getPartTaxa(int part, MTree *tree, MTree *part_tree, NodeVector &part_taxa){
    
    NodeVector taxa_nodes;
    tree->getTaxa(taxa_nodes);
    
    part_taxa.resize(taxa_num);
    
    if(!flag_reorderAccordingToTree){
        reorderAccordingToTree(taxa_nodes);
    }
    
    Node *node;
    for(NodeVector::iterator it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        //cout<<(*it)->name<<" id = "<<(*it)->id<<endl;
        if(pr_ab_matrix[(*it)->id][part] == 1){
            node = part_tree->findLeafName((*it)->name);
            //cout<<"PREPARING PART_TAXA: part = "<<part<<"|leaf_id = "<<(*it)->id<<"|leaf_name = "<<(*it)->name<<endl;
            assert(node && "ERROR: The leaf is not found on partition tree!");
            part_taxa[(*it)->id]=node;
        }
    }
    
    bool check = false;
    if(check){
        //cout<<"GET_taxa_nodes_info...."<<endl;
        cout<<"Presence-absence info:"<<endl;
        for(int i=0; i<taxa_num; i++){
            cout<<pr_ab_matrix[i][part]<<" ";
        }
        cout<<endl;
        cout<<"Partition taxa info:"<<endl;
        for(int i=0; i<taxa_num; i++){
            if(part_taxa[i]){
                cout<<part_taxa[i]->name<<"("<<part_taxa[i]->id<<") ";
            }else{
                cout<<"NA ";
            }
        }
        cout<<endl;
    }
}

void PresenceAbsenceMatrix::reorderAccordingToTree(NodeVector taxa_nodes){

    //cout<<"BEFORE reordering according to the tree:"<<endl;
    //print_pr_ab_matrix();
    
    int i=0;
    int id=-1;
    
    vector<IntVector> aux_matrix;
    vector<string> aux_names;
    
    aux_matrix.resize(taxa_num);
    aux_names.resize(taxa_num);
    
    for(i=0; i<taxa_nodes.size();i++){
        id=findTaxonID(taxa_nodes[i]->name);
        //cout<<taxa_nodes[i]->name<<" id = "<<id<<endl;
        aux_matrix[taxa_nodes[i]->id]=pr_ab_matrix[id];
        aux_names[taxa_nodes[i]->id]=taxa_names[id];
    }
    
    pr_ab_matrix.clear();
    taxa_names.clear();
    
    for(i=0; i<taxa_num; i++){
        pr_ab_matrix.push_back(aux_matrix[i]);
        taxa_names.push_back(aux_names[i]);
    }
    
    //cout<<"AFTER reordering according to the tree:"<<endl;
    //print_pr_ab_matrix();
}

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree){
    
    //int id;
    vector<IntVector> sub_matrix;
    NodeVector taxa_nodes;
    NodeVector::iterator it;
    string taxon_name;
    tree->getTaxa(taxa_nodes);
    for(it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        taxon_name=(*it)->name;
        //id=getTaxonID_in_pr_ab_m(taxon_name);
    }
    
    return sub_matrix;
};

void PresenceAbsenceMatrix::getSubPrAbMatrix(vector<string> taxa_names_subset, PresenceAbsenceMatrix *submatrix, IntVector *parts){
 
    vector<string> not_found_taxon_names;
    bool found = false;
    
    int i,j,h;
    for(i=0; i<taxa_names_subset.size(); i++){
        //cout<<i<<": subset_taxon = "<<taxa_names_subset[i]<<endl;
        found = false;
        for(j=0; j<taxa_names.size(); j++){
            //cout<<j<<": taxon_name = "<<taxa_names[j]<<endl;
            if(taxa_names_subset[i]==taxa_names[j]){
                //cout<<"MATCH"<<endl;
                if(parts){
                    IntVector taxon_coverage;
                    for(IntVector::iterator k=parts->begin(); k<parts->end(); k++){
                        h = (*k);
                        taxon_coverage.push_back(pr_ab_matrix[j][h]);
                    }
                    submatrix->pr_ab_matrix.push_back(taxon_coverage);
                }else{
                    submatrix->pr_ab_matrix.push_back(pr_ab_matrix[j]);
                }
                submatrix->taxa_names.push_back(taxa_names[j]);
                found = true;
                break;
            }
        }
        if(!found){
            cout<<"Taxon "<<taxa_names_subset[i]<<" is not found in the presence-absence matrix..."<<endl;
            not_found_taxon_names.push_back(taxa_names_subset[i]);
        }
    }
    
    if(not_found_taxon_names.size()<taxa_names_subset.size()){
        submatrix->taxa_num = submatrix->taxa_names.size();
        submatrix->part_num = submatrix->pr_ab_matrix[0].size();
        cout<<"INFO: original matrix."<<endl;
        print_pr_ab_matrix();
        cout<<endl<<"INFO: a submatrix for "<<submatrix->taxa_num<<" taxa was extracted."<<endl;
        submatrix->print_pr_ab_matrix();
    }
}

void PresenceAbsenceMatrix::getSubPrAbMatrix(NodeVector taxon_nodes, PresenceAbsenceMatrix *submatrix, IntVector *parts){

    vector<string> taxon_names;
    for(NodeVector::iterator it = taxon_nodes.begin(); it<taxon_nodes.end(); it++){
        taxon_names.push_back((*it)->name);
    }
    
    getSubPrAbMatrix(taxon_names,submatrix,parts);
}
