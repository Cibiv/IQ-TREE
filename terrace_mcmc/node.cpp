/*
 * node.cpp
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#include <boost/range/adaptor/map.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/functional/hash.hpp>

#include <iostream>
#include <algorithm>
#include <map>
#include "node.h"
#include "parser.h"
#include "supermatrix.h"

using namespace std;

MCNode::MCNode()
{
}

MCNode* MCNode::addChild()
{
    MCNode* child = new MCNode();
    child->parent = this;

    if (left_child != NULL) {
        if (right_child != NULL) {
            if (isRoot() || true) {
                MCNode* extraNode = new MCNode();
                extraNode->parent = this;
                extraNode->left_child = right_child;
                extraNode->right_child = child;
                child->parent = extraNode; 
                right_child = extraNode;
            } else {
                cout << "polytomy " << name << endl;
                cout << " left: " << left_child->name << endl;
                cout << " right: " << right_child->name << endl;
                throw "";
            }            
        } else {
            right_child = child;
        }        
    } else {
        left_child = child;
    }

    return child;
}

bitset<MAX_LOCI>* MCNode::recalculateLoci()
{
    return NULL;

    /*if (loci->any()) {
        //return loci;
    }

    bitset<MAX_LOCI>* child_loci;

    if (left_child != NULL) {
        child_loci = left_child->recalculateLoci();
        loci->operator|=(*child_loci);
    }

    if (right_child != NULL) {
        child_loci = right_child->recalculateLoci();
        loci->operator|=(*child_loci);  
    }

    return loci;*/
}

bitset<MAX_TAXA>* MCNode::recalculateLeaves()
{
    if (!leaves_reset) {
        return leaves;
    }

    leaves_reset = false;
    leaves->reset();

    bitset<MAX_TAXA>* child_leaves;

    if (left_child != NULL) {        
        child_leaves = left_child->recalculateLeaves();
        leaves->operator|=(*child_leaves);
    }

    if (right_child != NULL) {
        child_leaves = right_child->recalculateLeaves();
        leaves->operator|=(*child_leaves);  
    }

    return leaves;
}

void MCNode::resetLeaves()
{
    if (isLeaf()) {
        return;
    }

    leaves_reset = true;

    left_child->resetLeaves();
    right_child->resetLeaves();
}

bitset<MAX_TAXA>* MCNode::getLeaves()
{
    recalculateLeaves();
    return leaves;    
}

bitset<MAX_LOCI>* MCNode::getLoci()
{
    /*recalculateLoci();*/
    //return loci;
    return NULL;
}

bool MCNode::isRoot()
{
    return is_root;//name == "__root__";
}

bool MCNode::isLeaf()
{
    return is_leaf; //name != "" && !isRoot();
}

void MCNode::setName(string new_name)
{
    name = new_name;

    is_root = false;
    is_leaf = false;

    if (name == "__root__") {
        is_root = true;
    } else if (name != "") {
        is_leaf = true;
    }
}

string MCNode::getName()
{
    return name.length() > 0 ? name : "(None)";
}

void MCNode::swapParents(MCNode* a, MCNode* b)
{
    MCNode* a_parent = a->parent;
    MCNode* b_parent = b->parent;

    a->parent = b_parent;
    b->parent = a_parent;

    if (a_parent->left_child == a) {
        a_parent->left_child = b;
    } else {
        a_parent->right_child = b;
    }

    if (b_parent->left_child == b) {
        b_parent->left_child = a;
    } else {
        b_parent->right_child = a;
    }
}

void MCNode::performNNI(MCNode* nni[4], bool variant, bool noreset)
{
    //cout << "A: " << nni[0]->toString() << endl;
    //cout << "B: " << nni[1]->toString() << endl;
    //cout << "C: " << nni[2]->toString() << endl;
    //cout << "D: " << nni[3]->toString() << endl;
    //cout << variant << endl;
    //cout << print() << endl;
    //cout << "--" << endl;

    if (variant == true) {
        swapParents(nni[0], nni[2]);
    } else {
        swapParents(nni[1], nni[2]);
    }

    if (!noreset) {
        resetNNICache();
        resetLeaves(); 
    }    
}

string MCNode::getFirstLeaf()
{
    if (isLeaf()) {
        return name;
    }

    return left_child->getFirstLeaf();
}

MCNode* MCNode::clone()
{
    MCNode* copy = new MCNode(*this);

    copy->leaves = new bitset<MAX_TAXA>(*this->getLeaves());
    copy->resetNNICache();

    copy->num_nnis = -1;

    if (left_child != NULL) {
        copy->left_child = left_child->clone();
        copy->left_child->parent = copy;
    }

    if (right_child != NULL) {
        copy->right_child = right_child->clone();
        copy->right_child->parent = copy;
    }

    return copy;
}

string MCNode::toString()
{
    if (isRoot()) {
        return "root";
    }

    string out = "MCNode, ";

    if (isLeaf()) {
        out += "label ";
        out += name;
        out += ", ";
    }

    out += "leaves: ";    

    out += getLeaves()->to_string();

    return out;
}

bool MCNode::hasCommonLeaves(bitset<MAX_TAXA>* a, bool complement, MCNode* n)
{     
    bitset<MAX_TAXA> _leaves;

    if (complement) {
        _leaves = bitset<MAX_TAXA>(*n->getLeaves());        
        _leaves.flip();                
    } else {
        _leaves = bitset<MAX_TAXA>(*getLeaves());
    }

    _leaves &= *a;

    return _leaves.any();    
}

void MCNode::printNeighbours()
{
    string nwk_parent = printSorted();
    getAllNNIs();

    for (int i=0; i<num_nnis; i++) {        
        MCNode* nni[4];            

        for (int j=0; j<4; j++) {
            nni[j] = nni_cache[i*4 + j];    
        }

        for (int j=0; j<2; j++) {
            performNNI(nni, j == 0 ? true : false, true);
            cout << printSorted() << endl;
            performNNI(nni, j == 0 ? true : false, true);
        }        
    }   
}

void MCNode::_getAllNNIs(int &result_count, MCNode* result[4*(MAX_TAXA-3)])
{
    recalculateLeaves();
    //recalculateLoci();

    if (isLeaf()) {
        return;
    }

    bool nniComplementD = false;

    if (isRoot()) {
        if (right_child->isLeaf()) {
            MCNode* tmp = right_child;
            right_child = left_child;
            left_child = tmp;
        }

        left_child->_getAllNNIs(result_count, result);
        right_child->_getAllNNIs(result_count, result);        
        return;
    }

    if (parent->isRoot()) {
        if (this == parent->left_child) {            
            result[result_count*4] = left_child;
            result[result_count*4+1] = right_child; 
            result[result_count*4+2] = parent->otherChild(this)->left_child;
            result[result_count*4+3] = parent->otherChild(this)->right_child;                                
        } else {
            left_child->_getAllNNIs(result_count, result);
            right_child->_getAllNNIs(result_count, result); 
            return;           
        }
    } else if (parent->parent->isRoot()) {
        result[result_count*4] = left_child;
        result[result_count*4+1] = right_child; 
        result[result_count*4+2] = parent->otherChild(this);
        result[result_count*4+3] = parent->parent->otherChild(parent);     
    } else {
        result[result_count*4] = left_child;
        result[result_count*4+1] = right_child; 
        result[result_count*4+2] = parent->otherChild(this); 
        result[result_count*4+3] = parent->parent;

        nniComplementD = true;
    } 

    result_count++;

    left_child->_getAllNNIs(result_count, result);
    right_child->_getAllNNIs(result_count, result);
}

void MCNode::getAllNNIs()
{
    if (num_nnis == -1) {
        num_nnis = 0;
        _getAllNNIs(num_nnis, nni_cache);
    }    
}

{void MCNode::_getOnTerraceNNIs(Supermatrix* supermatrix, int &result_count, MCNode* result[4*(MAX_TAXA-3)])

    recalculateLeaves();
    //recalculateLoci();

    if (isLeaf()) {
        return;
    }

    bool nniComplementD = false;

    if (isRoot()) {
        if (right_child->isLeaf()) {
            MCNode* tmp = right_child;
            right_child = left_child;
            left_child = tmp;
        }

        left_child->_getOnTerraceNNIs(supermatrix, result_count, result);
        right_child->_getOnTerraceNNIs(supermatrix, result_count, result);        
        return;
    }

    if (parent->isRoot()) {
        if (this == parent->left_child) {            
            result[result_count*4] = left_child;
            result[result_count*4+1] = right_child; 
            result[result_count*4+2] = parent->otherChild(this)->left_child;
            result[result_count*4+3] = parent->otherChild(this)->right_child;                                
        } else {
            left_child->_getOnTerraceNNIs(supermatrix, result_count, result);
            right_child->_getOnTerraceNNIs(supermatrix, result_count, result); 
            return;           
        }
    } else if (parent->parent->isRoot()) {
        result[result_count*4] = left_child;
        result[result_count*4+1] = right_child; 
        result[result_count*4+2] = parent->otherChild(this);
        result[result_count*4+3] = parent->parent->otherChild(parent);     
    } else {
        result[result_count*4] = left_child;
        result[result_count*4+1] = right_child; 
        result[result_count*4+2] = parent->otherChild(this); 
        result[result_count*4+3] = parent->parent;

        nniComplementD = true;
    } 

    static bitset<MAX_TAXA>* bitsets[MAX_LOCI];
    static int bitset_count = 0;
    bitset_count    = 0;

    for(auto&& i : supermatrix->taxon_by_locus | boost::adaptors::map_values){
        bitsets[bitset_count] = i;
        bitset_count++;
    }

    bool is_on_terrace = true;
 
    for (int i=0; i<bitset_count; i++) {
        if (
            result[result_count*4]->hasCommonLeaves(bitsets[i]) &&
            result[result_count*4+1]->hasCommonLeaves(bitsets[i]) &&
            result[result_count*4+2]->hasCommonLeaves(bitsets[i])  &&
            result[result_count*4+3]->hasCommonLeaves(bitsets[i], nniComplementD, parent)
        ) {
            is_on_terrace = false;
            break;              
        }
    }

    if (is_on_terrace) {         
        result_count++;
    }            

    left_child->_getOnTerraceNNIs(supermatrix, result_count, result);
    right_child->_getOnTerraceNNIs(supermatrix, result_count, result);
}

void MCNode::getOnTerraceNNIs(Supermatrix* supermatrix)
{
    if (num_nnis == -1) {
        num_nnis = 0;
        _getOnTerraceNNIs(supermatrix, num_nnis, nni_cache);
    }    
}

void MCNode::replaceChild(MCNode* child, MCNode* replacement)
{
    if (parent->isRoot()) {
        replacement = parent->otherChild(this);        
    }  else {
        parent->replaceChild(this, parent->parent);
    } 

    if (left_child == child) {
        left_child = replacement;
        replacement->parent = this;
    } else if (right_child == child) {
        right_child = replacement;
        replacement->parent = this;
    } else {
        throw "invalid";
    }    
}

void MCNode::setRoot(bool root)
{
    is_root = root;
}

void MCNode::root(MCNode* root, string taxon)
{
    if (isLeaf()) {
        if (name == taxon and !parent->isRoot()) {
            MCNode* p = this->parent;
            MCNode* previous = this;

            p->replaceChild(previous, p->parent);

            root->left_child = this;
            root->right_child = this->parent; 

            this->parent->parent = root;
            this->parent = root;        
        }        

        return;
    }

    left_child->root(root, taxon);

    if (right_child != NULL) {
        right_child->root(root, taxon);
    }

    resetLeaves();
    resetNNICache();
}

void MCNode::resetNNICache()
{
    if (num_nnis != -1) {
        num_nnis = -1;
    }

    if (left_child != NULL) {
        left_child->resetNNICache();
    }

    if (right_child != NULL) {
        right_child->resetNNICache();
    }      
}

void MCNode::sort()
{
    if (isLeaf()) {
        return;
    }

    MCNode* lc = left_child;
    MCNode* rc = right_child;

    lc->sort();
    rc->sort();

    bitset<MAX_TAXA>* l = lc->getLeaves();
    bitset<MAX_TAXA>* r = rc->getLeaves();

    for (int i=0; i<MAX_TAXA; i++) {        
        if (l->test(i) != r->test(i)) {
            if (r->test(i)) {
                right_child = lc;
                left_child = rc;
            }

            break;
        }        
    }
}

void MCNode::setTaxonLoci(map<int, bitset<MAX_LOCI>*>* locus_by_taxon)
{
  /*  if (isLeaf()) {
        loci = locus_by_taxon->at(name_index);
    } else {
        left_child->setTaxonLoci(locus_by_taxon);
        right_child->setTaxonLoci(locus_by_taxon);
    } */   
}

void MCNode::getNNI(int no, MCNode* nni[4])
{
    nni[0] = nni_cache[no];
    nni[1] = nni_cache[no+1];
    nni[2] = nni_cache[no+2];
    nni[3] = nni_cache[no+3];    
}

void MCNode::setTaxonIndices(map<std::string, int>* index)
{
    if (isLeaf()) {
        //cout << "--" << endl;
        //cout << name << endl;
        name_index = index->at(name);
        //cout << name_index << endl;
        //cout << "--" << endl;
        leaves->set(name_index);
        leaves_reset = false;
    } else {
        left_child->setTaxonIndices(index);
        right_child->setTaxonIndices(index);
    }
}

MCNode* MCNode::otherChild(MCNode* child)
{
    if (child == left_child) {
        return right_child;
    } 

    return left_child;    
}

MCNode::~MCNode()
{
    delete leaves;

    //if (!isLeaf()) {
    if (left_child != NULL)
        delete left_child;

    if (right_child != NULL)
        delete right_child;
}

void MCNode::deleteChild(MCNode* child) 
{
    MCNode* other = otherChild(child);

    if (left_child == child) {
        left_child = right_child;        
    }

    right_child = NULL;

    delete child;    

    if (!isRoot()) {
        left_child = NULL;
        right_child = NULL;

        if (other == NULL) { // both children null
           // cout << "deleting child" << endl;
            parent->deleteChild(this);
        } else {
        //    cout << "replacing child" << endl;
            parent->replaceChild2(this, other);
        }    
    }
}

void MCNode::replaceChild2(MCNode* child, MCNode* replacement)
{
    replacement->parent = this;

    if (left_child == child) {
        left_child = replacement;
    } else {
        right_child = replacement;
    }

    delete child;
}

void MCNode::getLeafNodes(vector<MCNode*> &nodes)
{
    if (left_child->isLeaf()) {
        nodes.push_back(left_child);
    } else {
        left_child->getLeafNodes(nodes);
    }

    if (right_child->isLeaf()) {
        nodes.push_back(right_child);
    } else {
        right_child->getLeafNodes(nodes);
    }    
}

void MCNode::getSubtree(set<string> &taxa)
{
    vector<MCNode*> nodes;
    getLeafNodes(nodes);

    for (std::vector<MCNode*>::iterator it = nodes.begin() ; it != nodes.end(); ++it) {
        cout << (*it)->name << endl;
        if (taxa.find((*it)->name) == taxa.end()) {
            (*it)->parent->deleteChild(*it);
        } 
    }

  //  cout << "A" << endl;

    if (left_child != NULL && right_child != NULL) {
        return;
    }

    MCNode* onlyChild = left_child;
    if  (onlyChild == NULL) {
        onlyChild = right_child;
    }

    if (onlyChild->isLeaf()) {
        return;
    }

    left_child = NULL;
    right_child = NULL;

    left_child = onlyChild->left_child;
    left_child->parent = this;

    if (onlyChild->right_child != NULL) {
        right_child = onlyChild->right_child;
        right_child->parent = this;
    }

    onlyChild->left_child = NULL;
    onlyChild->right_child = NULL;

    delete onlyChild;
}

void MCNode::getTaxonNames(set<string> &names)
{
    if (isLeaf()) {
        names.insert(name);
    } else {
        left_child->getTaxonNames(names);
        right_child->getTaxonNames(names);
    }
}

string MCNode::getSupertreeString(Supermatrix* sm)
{
    string output = "";

    for (int i=0; i<sm->total_loci; i++) { 
        MCNode* copy = clone();
        set<string> taxa = sm->getPartitionTaxa(i);

        copy->getSubtree(taxa);
        string root_label = *min_element(taxa.begin(), taxa.end());

        copy->root(this, root_label);

        output += copy->printSorted();
        output += ";";

        delete copy;
    }

    return output;
}

std::string get_sha1(const std::string& p_arg)
{

    boost::uuids::detail::sha1 sha1;
    sha1.process_bytes(p_arg.data(), p_arg.size());
    unsigned hash[5] = {0};
    sha1.get_digest(hash);

    // Back to string
    char buf[41] = {0};

    for (int i = 0; i < 5; i++)
    {
        std::sprintf(buf + (i << 3), "%08x", hash[i]);
    }

    return std::string(buf);
}

string MCNode::getTerraceId(Supermatrix* sm)
{
    string str = getSupertreeString(sm);
    string hash = get_sha1(str);
    str += " ";
    str += hash;

    return str;
}

string MCNode::getSupertreeString(string supertreeString)
{
    string output = print(false);
    output += ";";

    int i = 0;
    string first = "";
    while (supertreeString[i] != ';') { first += supertreeString[i]; i++; }
    first += ";";

    cout << supertreeString << endl;

    //cout << p->rmweights(first) << endl;

    int partindex = 0;
    while (i < supertreeString.length() - 1) {
        string nwk = "";

        do {
            i++;
            nwk += supertreeString[i];
        } while (supertreeString[i] != ';');

        //cout << "NWK SUPERSTRINGPART " << nwk << endl;

        Parser* p = new Parser();
        //cout << "NWK " << nwk << endl;
        MCNode* parsed = p->from_newick(nwk); // induced subtree
        //cout << "ROOTED " << p->root(nwk) << endl;
        set<string> taxa;
        parsed->getTaxonNames(taxa);

        MCNode* copy = clone(); // copy of comprehensive
        //cout << "A" << endl;
        copy->getSubtree(taxa); 
        //cout << "B" << endl;
        string cpp = copy->print(false);
        cpp += ";"; 

        //cout << "--" << p->rmweights(cpp) << endl;

        output += cpp;
        partindex++;

        delete parsed;
        delete copy;
        delete p;
    }

    return output;   
}

void MCNode::replaceLeafNames(std::map<int,int> aln)
{
    vector<MCNode*> nodes;
    getLeafNodes(nodes);

    for (std::vector<MCNode*>::iterator it = nodes.begin() ; it != nodes.end(); ++it) {
        int from = stoi((*it)->name);

        if (aln.find(from) == aln.end()) {
            cout << "NOT FOUND ID " << from << endl;
            throw "";
        }

        (*it)->name = to_string(aln[from]);
    }        
}

string MCNode::getSupertreeString(string supertreeString, std::map<int, std::map<int,int>> aln, std::map<int, std::map<int,int>> revaln)
{
    //cout << "SSTRING " << supertreeString << endl;
    string output = print(false);
    output += ";";

    int i = 0;
    string first = "";
    while (supertreeString[i] != ';') { first += supertreeString[i]; i++; }
    first += ";";

    //cout << p->rmweights(first) << endl;

    int partindex = 0;
    while (i < supertreeString.length() - 1) {
        string nwk = "";

        do {
            i++;
            nwk += supertreeString[i];
        } while (supertreeString[i] != ';');

        Parser* p = new Parser();
        //cout << "NWK " << nwk << endl;
        //cout << "ROOTED " << p->root(nwk) << endl;
        MCNode* parsed = p->from_newick(nwk); // induced subtree
        //cout << "DONE" << endl;

       // cout << "WITH COMP " << p->rmweights(parsed->print()) << endl;
        parsed->replaceLeafNames(aln[partindex]); // from partition ids to comprehensive ids
      //  cout << "WITH PART " << p->rmweights(parsed->print()) << endl;

        set<string> taxa;
        parsed->getTaxonNames(taxa);

        MCNode* copy = clone(); // copy of comprehensive
        //cout << "A" << endl;
        copy->getSubtree(taxa); 
        //cout << "B" << endl;

      //  cout << "WITH PART " << p->rmweights(copy->print()) << endl;
        copy->replaceLeafNames(revaln[partindex]);
      //  cout << "WITH COMP " << p->rmweights(copy->print()) << endl;
      //  cout << "--" << endl;

        string cpp = copy->print(false);
        cpp += ";";         

        output += cpp;
        partindex++;

        delete p;
        delete parsed;
        delete copy;
    }

    return output;
}

void MCNode::debug(string indent)
{
    string newindent = indent;
    newindent += "  ";

    if (name != "")
        cout << name << endl;
    else
        cout << "(internal)" << endl;

    if (parent == NULL) {
        cout << indent << "PARENT: NULL" << endl;
    } else {
        cout << indent << "PARENT: " << parent->toString() << endl;
    }

    if (left_child == NULL) {
        cout << indent << "Left: NULL" << endl;
    } else {
        cout << indent << "Left: ";
        left_child->debug(newindent);
    }

    if (right_child == NULL) {
        cout << indent << "Right: NULL" << endl;
    } else {
        cout << indent << "Right: ";
        right_child->debug(newindent);
    }
}

string MCNode::printSorted(bool rooted, bool weights)
{    
    sort();
    return print(rooted, weights);     
}

string MCNode::print(bool rooted, bool weights)
{
    weights = true;
    
    if (isLeaf()) {
        string out = name;

        //if (weights)
           // out += ":0.0000010883";

        return out;
    }

    string out = "(";
    string left, right;

    if (left_child != NULL) {
        if (right_child != NULL && !rooted && right_child->isLeaf() && !left_child->isLeaf()) {
            left = left_child->left_child->print(true, weights);
            left += ",";
            left += left_child->right_child->print(true, weights);            
        } else {
            left = left_child->print(true, weights);
        }    
    } else {
        left = "";
    }

    if (right_child != NULL) {
        if (isRoot() && !rooted && !right_child->isLeaf()) {
            right = right_child->left_child->print(true, weights);
            right += ",";
            right += right_child->right_child->print(true, weights);
        } else {
            right = right_child->print(true, weights);
        }
    } else {
        right = "";
    }

    if (left != "" && right != "") {
        left += ",";
    }

    if (left != "")
        out += left;

    if (right != "")
        out += right;

    out += ')';

    if (!isRoot() && weights) {
        //out += ":0.0000010883";   
    }

    if (isRoot()) {
        out += ";";
    }

    return out;
}
