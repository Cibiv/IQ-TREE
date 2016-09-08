//
//  cladeanalysis.cpp
//  iqtree
//
//  Created by Olga on 31/08/16.
//
//

#include "cladeAnalysis.h"

/**
 *  Here we are interested in analysis of clades. The input files are the tree and the list of taxa. The question is whether the taxa from the list form a separate clade or whether they are mixed with some other taxa, i.e. if a common anscestor is also an anscestor of some other taxa on the tree. As an output it would be nice to have yes/no and if no, the list of taxa that belong to the same clade.
 **/

CladeAnalysis::CladeAnalysis(IQTree* tree) : MTree(){
    initCladeAnalysis(tree);
};

CladeAnalysis::~CladeAnalysis(){
};

void CladeAnalysis::initCladeAnalysis(IQTree* tree){
    
    minCladeSize = tree->leafNum;
    
    // Reading taxa to search the clade for
    readInputTaxa(Params::getInstance().clade_analysis_infile);
    
    // Check whether all the taxa from the list are present on the tree
    for(int i = 0; i < taxaName.size(); i++){
        if(!tree->findLeafName(taxaName[i])){
            cout<<"ERROR: Taxon "<<taxaName[i]<<" is not present on the tree!!"<<endl;
            cout<<"Check the spelling. Stopping this analysis."<<endl;
            exit(0);
        } else {
            taxaNameID.push_back(tree->findLeafName(taxaName[i])->id);
        }
    }
    
};

void CladeAnalysis::startCladeAnalysis(IQTree* tree){
    
    int i, j;
    
    cout<<"========================================================="<<endl;
    cout<<"         Starting Clade analysis.."<<endl;
    cout<<"========================================================="<<endl;
    
    int taxaNameNUM = taxaName.size();
    
    NodeVector branch1, branch2;
    tree->getBranches(branch1, branch2);
    
    // A loop over all A|B present on tree T
    vector<int> taxaA, taxaB;
    int taxaAsize = 0, taxaBsize = 0;
    
    bool foundALL = FALSE;
    bool foundSOME = FALSE;
    
    string taxa_set;
    for(i = 0; i < tree->leafNum; i++)
        taxa_set.push_back(0);
    
    for(i = 0; i != branch1.size(); i++){
        tree->getTaxaID(taxaA,branch1[i],branch2[i]);
        tree->getTaxaID(taxaB,branch2[i],branch1[i]);
        
        taxaAsize = taxaA.size();
        taxaBsize = taxaB.size();
        
        // We consider only non-trivial splits
        if(taxaAsize > 1 && taxaBsize > 1){
            
            if(taxaAsize <= taxaBsize){
                checkClade(taxaA, &foundALL, &foundSOME);
            } else {
                checkClade(taxaB, &foundALL, &foundSOME);
            }
            
            if(taxaAsize >= taxaNameNUM){
                if(taxaAsize < minCladeSize){
                                        if(foundALL){
                        minCladeSize = taxaAsize;
                        minCladeSpecies = taxaAname;
                        for(i = 0; i < taxaAsize; i++)
                            taxa_set[taxaA[i]]=1;
                        ((MTree*)tree)->copyTree(&minCladeSubtree, taxa_set);
                        for(i = 0; i < taxaAsize; i++)
                            taxa_set[taxaA[i]]=0;
                    }
                    
                }
            }
            
            if(!foundALL && taxaBsize >= taxaNameNUM && taxaBsize < minCladeSize){
                i = 0;
                foundALL = TRUE;
                while(foundALL && i<taxaNameNUM){
                    foundALL = FALSE;
                    for(j = 0; j < taxaBsize; j++){
                        if(strcmp(taxaName[i].c_str(),taxaBname[j].c_str()) == 0){
                            foundALL = TRUE;
                            cout<<"Found species "<<i<<": "<<taxaName[i]<<endl;
                            break;
                        }
                    }
                    if(!foundALL)
                        cout<<" ---> Did not find species "<<i<<": "<<taxaName[i]<<endl;
                    i++;
                }
                if(foundALL){
                    minCladeSize = taxaBsize;
                    minCladeSpecies = taxaBname;
                    for(i = 0; i < taxaBsize; i++)
                        taxa_set[taxaB[i]]=1;
                    ((MTree*)tree)->copyTree(&minCladeSubtree, taxa_set);
                    for(i = 0; i < taxaBsize; i++)
                        taxa_set[taxaB[i]]=0;
                }
                
                
            }

            
        }// END of We consider only non-trivial splits
    }
    
    
}

void CladeAnalysis::readInputTaxa(const char* infile){
    ifstream in;
    cout<<endl<<"-----------------------------------------------------"<<endl;
    cout<<"Reading taxa of interest from "<<infile<<endl;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        readInputTaxa(in);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
}

void CladeAnalysis::readInputTaxa(istream &in) {
    int i;
    string str,name;
    while (getline(in, str)) {
        stringstream ss(str);
        getline(ss,name,':');
        taxaName.push_back(name);
    }
    for(i=0; i<taxaName.size();i++)
        cout<<taxaName.at(i)<<endl;
}

void CladeAnalysis::checkClade(vector<string> *taxaSplit, bool *foundALL, bool *foundSOME){
}

void CladeAnalysis::setMinClade(IQTree *tree, vector<string> *taxaSplit){
}
