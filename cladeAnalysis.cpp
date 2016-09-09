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
    
    taxaNameNUM = taxaName.size();
    // Check whether all the taxa from the list are present on the tree
    for(int i = 0; i < taxaNameNUM; i++){
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
    
    int i;
    
    cout<<"========================================================="<<endl;
    cout<<"         Starting Clade analysis.."<<endl;
    cout<<"========================================================="<<endl;
    
    NodeVector branch1, branch2;
    tree->getBranches(branch1, branch2);
    
    // A loop over all A|B present on tree T
    vector<int> taxaA, taxaB;
    int taxaAsize = 0, taxaBsize = 0;
    
    bool foundALL = false;
    bool foundSOME = false;
    
    for(i = 0; i != branch1.size(); i++){
        tree->getTaxaID(taxaA,branch1[i],branch2[i]);
        tree->getTaxaID(taxaB,branch2[i],branch1[i]);
        
        taxaAsize = taxaA.size();
        taxaBsize = taxaB.size();
        
        // We consider only non-trivial splits
        if(taxaAsize > 1 && taxaBsize > 1){
            if(taxaAsize <= taxaBsize){
                checkClade(&taxaA, &foundALL, &foundSOME);
                if(foundALL && taxaAsize < minCladeSize){
                    setMinClade(tree, &taxaA);
                } else if (!foundSOME && taxaBsize < minCladeSize) {
                    setMinClade(tree, &taxaB);
                }
            } else {
                checkClade(&taxaB, &foundALL, &foundSOME);
                if(foundALL && taxaBsize < minCladeSize){
                    setMinClade(tree, &taxaB);
                } else if (!foundSOME && taxaAsize < minCladeSize) {
                    setMinClade(tree, &taxaA);
                }
            }
        }
    }
    
    printResultsCA();
    exit(0);
    
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

void CladeAnalysis::checkClade(vector<int> *taxaSplit, bool *foundALL, bool *foundSOME){
    int i = 0, j = 0;
    int taxaSplitSize = taxaSplit->size();
    
    *foundALL = true;
    
    //  We continue checking all species in taxaSplit till either:
    //  (i)     so far we found all taxa: foundALL = true, or
    //  (ii)    so far all taxa are not present in the split (foundALL = false, so far), but we don't know if some other is still present (foundSOME = false, so far)
    
    while(i < taxaNameNUM && (*foundALL || (!(*foundALL) && !(*foundSOME)))){
        *foundALL = false;
        for(j = 0; j < taxaSplitSize; j++){
            if(taxaNameID[i] == taxaSplit->at(j)){
                *foundALL = true;
                *foundSOME = true;
                cout<<"Found species "<<i<<": name = "<<taxaName[i]<<"; id = "<<taxaNameID[i]<<endl;
                break;
            }
        }
        if(!(*foundALL))
            cout<<" ---> Did not find species "<<i<<": name = "<<taxaName[i]<<"; id = "<<taxaNameID[i]<<endl;
        
        // For control only
        cout<<"Species "<<i<<": foundALL = "<<*foundALL<<"; foundSOME"<<*foundSOME<<endl;
        
        i++;
    }

}

void CladeAnalysis::setMinClade(IQTree *tree, vector<int> *taxaSplit){
    int i;
    
    assert(taxaSplit->size()<minCladeSize);
    
    cout<<"Updating minClade...-------------------------------"<<endl;
    cout<<"     old size of minClade = "<<minCladeSize<<endl;
    
    minCladeSize = taxaSplit->size();
    
    cout<<"     new size of minClade = "<<minCladeSize<<endl;
    cout<<"---------------------------------------------------"<<endl;
    
    for(i = 0; i < minCladeSize; i++){
        minCladeSpecies.push_back(tree->findNodeID(taxaSplit->at(i))->name);
    }
    string taxa_set;
    for(i = 0; i < tree->leafNum; i++)
        taxa_set.push_back(0);
    for(i = 0; i < minCladeSize; i++)
        taxa_set[taxaSplit->at(i)]=1;
    ((MTree*)tree)->copyTree(this, taxa_set);
}

void CladeAnalysis::printResultsCA(){
    
    bool found = false;
    int i, j = 0;
    string out_file = Params::getInstance().out_prefix;
    out_file += ".ca.details";
    
    ofstream out;
    out.exceptions(ios::failbit | ios::badbit);
    out.open((char*)out_file.c_str(),std::ofstream::out);
    out<<"minCladeSize "<<minCladeSize<<" minPossible "<<taxaNameNUM<<" ? "<<minPossible()<<endl;
    out.close();
    
    out_file = Params::getInstance().out_prefix;
    out_file += ".ca.minClade.subtree";
    this->printTree((char*)out_file.c_str());
    
    if(!minPossible()){
        out_file = Params::getInstance().out_prefix;
        out_file += ".ca.minClade.species";
        out.open((char*)out_file.c_str(),std::ofstream::out);
    
        for(i = 0; i < minCladeSpecies.size(); i++){
            found = false;
            
            while(!found && j < taxaNameNUM){
                if(this->findLeafName(minCladeSpecies.at(i))->id == taxaNameID.at(j))
                    found = true;
                j++;
            
            }
            if(!found)
                out<<minCladeSpecies.at(i)<<endl;
        }
        out.close();
    }
    
}

bool CladeAnalysis::minPossible(){
    if(minCladeSize == taxaNameNUM){
        return true;
    } else {
        return false;
    }
}

















