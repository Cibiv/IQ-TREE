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
#ifndef TERRACE_H2 //why _H2??
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
	 *  @param tree tree    //why tree and not the root?
	 *  @param saln superalignment
	 */
    Terrace(MCNode* root, Supermatrix* sm);

    /**
     *  @return The terrace size
     */
    uint64_t getSize();

    uint64_t size;
   
    std::vector<std::string> getSample(int size);

    bool isTrivial();

	void init();

    void printTrees();

    ~Terrace();

    MCNode* root;

    Supermatrix* sm;

    std::set<std::string> trees;

    void completeWalk();

    int burnin = -1;

    bool print_degrees = false;

    bool print_trees = false;

    bool print_unrooted = false;

private:

};


#endif
