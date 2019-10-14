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
#ifndef TERRACE_H2
#define TERRACE_H2

#include "supermatrix.h"
#include "srw.h"
#include "mhrw.h"
#include "parser.h"

/**
    A phylogenetic terrace
    @author Lukasz Reszczynski <lukasz.reszczynski@univie.ac.at>
*/
class Terrace {

public:
	/**
	 *  Constructor
	 *  @param tree tree
	 *  @param saln superalignment
	 */
    Terrace(MCNode* root, Supermatrix* sm);
    
    /*
     *  Destructor
     */
    ~Terrace();
    
    /*************************************************************************/

    /*
     * starting/representative tree on the terrace?
     */
    MCNode* root;
    
    /*
     * Presence-abscence matrix
     */
    Supermatrix* sm;
    
    /*
     *  Terrace size
     */
    uint64_t size;
    
    /*
     *  The set of trees on a terrace
     */
    std::set<std::string> trees;
    
    /*
     * burnin step for the MCMC sampling
     */
    int burnin = -1;
    
    /*
     * printing parameters
     */
    
    bool print_degrees = false;
    
    bool print_trees = false;
    
    bool print_unrooted = false;
    
    /*************************************************************************/
    
    void init();
    
    /**
     *  @return the terrace size
     */
    uint64_t getSize();
    
    /*
     *  get a sample from the terrace
     */
   
    std::vector<std::string> getSample(int size);

    /*
     *  check if the terrace is a trivial one, contains only one tree
     */
    bool isTrivial();
    
    /*
     *  Print out a sample of trees from the terrace
     */

    void printTrees();

    /*
     *  Perform a walk on the terrace
     */
    void completeWalk();
    

private:

};


#endif
