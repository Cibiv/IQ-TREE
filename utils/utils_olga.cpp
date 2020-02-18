//
//  utils_olga.cpp
//  iqtree
//
//  Created by Olga on 18.02.20.
//

#include "utils_olga.hpp"

void printSplitLength(Params &params){
    
    MTree tree(params.user_file, params.is_rooted);
    
    string out_split_len = params.out_prefix;
    out_split_len += ".split_len";
    ofstream out;
    out.open(out_split_len.c_str());
    
    NodeVector branch1, branch2;
    tree.getBranches(branch1, branch2);
    
    
    int i=0, j=0;
    
    // A loop over all A|B present on tree T
    for(i=0; i != branch1.size(); i++){
        vector<string> taxaAname, taxaBname;
        string splitA, splitB;
        
        tree.getTaxaName(taxaAname,branch1[i],branch2[i]);
        tree.getTaxaName(taxaBname,branch2[i],branch1[i]);
        
        //splitA.clear();
        for(j=0; j<taxaAname.size(); j++){
            splitA+=" "+taxaAname[j];
        }
        
        //splitB.clear();
        for(j=0; j<taxaBname.size(); j++){
            splitB+=" "+taxaBname[j];
        }
        
        if(taxaAname.size()<=taxaBname.size()){
            cout<<splitA<<"|"<<splitB<<"|"<<(branch1[i]->findNeighbor(branch2[i]))->length<<endl;
            out<<splitA<<"|"<<splitB<<"|"<<(branch1[i]->findNeighbor(branch2[i]))->length<<endl;
        }else{
            cout<<splitB<<"|"<<splitA<<"|"<<(branch1[i]->findNeighbor(branch2[i]))->length<<endl;
            out<<splitB<<"|"<<splitA<<"|"<<(branch1[i]->findNeighbor(branch2[i]))->length<<endl;
        }
        
    }
    
    out.close();
}
