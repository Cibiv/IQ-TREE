/*
 * utils.h
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#pragma once

#include <string>
#include <set>
#include <vector>

#include "parser.h"
#include "srw.h"
#include "mhrw.h"
#include "supermatrix.h"
#include "terrace.h"

class Utils {
    std::vector<std::string> getSample(std::string alignment_file, std::string tree_file, int size, int burnin);

    std::vector<std::string> getSample(Supermatrix* sm, std::string newick, int size, int burnin);

    std::vector<std::string> getSample(Supermatrix* sm,  MCNode* root, int size, int burnin);
};

