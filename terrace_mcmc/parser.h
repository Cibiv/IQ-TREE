/*
 * parser.h
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "node.h"

class Parser {
public:
    Parser();

    void load_node(MCNode* node);

    void load_children(MCNode* parent);

    bool in_range(char c, char from, char to);

    bool is_namechar(char c);

    std::string rmweights(std::string nwk);

    MCNode* from_file(std::string filename);

    MCNode* from_newick(std::string newick);

    std::string consume_name();

protected:
    const char* character;
};
