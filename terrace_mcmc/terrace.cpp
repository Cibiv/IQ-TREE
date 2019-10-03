/***************************************************************************
 *   Copyright (C) 2018 by Lukasz Reszczynski                              *
 *   lukasz.reszczynski@univie.ac.at                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "terrace.h"

using namespace std;

Terrace::Terrace(MCNode* root, Supermatrix* sm)
{
    root->setTaxonIndices(&sm->index);

    this->root = root;
    this->sm = sm;

    init();
}

void Terrace::init()
{

}

uint64_t Terrace::getSize()
{
    return size;
}

void Terrace::completeWalk()
{
    vector<string> queue;

    MCNode* tree_obj;
    string nwk;
    int num_nnis;

    string nwk_parent = root->printSorted();

    queue.push_back(nwk_parent);
    trees.insert(nwk_parent);

    Parser* p = new Parser();

    while (queue.size() > 0) {
        nwk_parent = queue.back();
        queue.pop_back();

        tree_obj = p->from_newick(nwk_parent);

        tree_obj->setTaxonIndices(&sm->index);
        tree_obj->getOnTerraceNNIs(sm);

        num_nnis = tree_obj->num_nnis;

        if (print_degrees) {
            cout << nwk_parent << " " << num_nnis << endl;
        } else if (print_trees) {
            cout << nwk_parent << endl;
        }

        for (int i=0; i<num_nnis; i++) {        
            MCNode* nni[4];            

            for (int j=0; j<4; j++) {
                nni[j] = tree_obj->nni_cache[i*4 + j];    
            }

            for (int j=0; j<2; j++) {
                tree_obj->performNNI(nni, j == 0 ? true : false);

                nwk = tree_obj->printSorted();

                if (trees.find(nwk) == trees.end()) {
                    queue.push_back(nwk);
                    trees.insert(nwk);
                }

                tree_obj->performNNI(nni, j == 0 ? true : false);
            }
        }

        delete tree_obj;
    }

    delete p;

    size = trees.size();
}

void Terrace::printTrees()
{
    for (std::set<string>::iterator it = trees.begin() ; it != trees.end(); ++it) {
        cout << *it << endl;
    }
}

bool Terrace::isTrivial()
{
    root->getOnTerraceNNIs(sm);
    return (root->num_nnis == 0);
}

vector<string> Terrace::getSample(int size)
{
    SRW* srw = new MHRW();
    srw->sm = sm;
    srw->burnin = (burnin > -1 ? burnin : atoi(getenv("BURNIN")));
    srw->sampling_freq = 1;
    srw->sample_size = 1;
    srw->command = "print_distinct";
    srw->print_degrees = print_degrees;
    srw->print_trees = print_trees;

    vector<string> sample;

    for (int i = 0; i < size; i++) { 
        if (srw->state) {
            delete srw->state;
        }   

        srw->state = root->clone();
        srw->state->resetNNICache();        

        string tree = srw->getSample()->at(0);
        sample.push_back(tree);
    }
    
    delete srw;

    return sample;
}

Terrace::~Terrace()
{
    //delete sm;
    delete root;    
}


