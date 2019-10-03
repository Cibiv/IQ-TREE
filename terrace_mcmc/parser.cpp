/*
 * parser.cpp
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */

#include "parser.h"

using namespace std;

Parser::Parser() 
{
    
}

bool Parser::in_range(char c, char from, char to)
{
    return c >= from && c <= to;
}

bool Parser::is_namechar(char c)
{        
    return in_range(c, '0', '9') ||
           in_range(c, 'A', 'Z') ||
           in_range(c, 'a', 'z') ||
           c == '_' ||
           c == '-' || 
           c == '.' || 
           c == '\'';
}

MCNode* Parser::from_file(string filename)
{
    ifstream f(filename.c_str());

    string newick;
    getline(f, newick);

    return from_newick(newick);
}

string Parser::rmweights(string in_newick)
{
    string newick = "";

    for (int i = 0; i<in_newick.length(); i++) {
        if (in_newick[i] == ':') {
            while (in_newick[i] != ')' && in_newick[i] != ',') i++;
        }        

        newick += in_newick[i];
    }

    return newick;
}

MCNode* Parser::from_newick(string newick)
{
    newick = rmweights(newick);

    character = newick.c_str();

    MCNode *root = new MCNode();
    
    load_node(root);
    root->setRoot();

    return root;    
}

string Parser::consume_name()
{
    std::string ret;

    while (is_namechar(*character))
    {
        if (*character != '\'') {
            ret += *character;
        }

        ++character;
    }

    return ret;
}

void Parser::load_node(MCNode* node)
{
    std::string name;
    //cout << "Loading node " << node->name << endl;

    switch (*character)
    {
        case '(':
            // We are nonleaf. Load new child.
            ++character;
            load_children(node); // leaves in a parent
            character++;
            name = consume_name();
            node->setName(name);
            break;
        case ',':
        case ')':
            //Allow nameless nodes: dont consume character.            
            //node->setName("");
            break;
        default:
            // We are leaf.
            name = consume_name();
            node->setName(name);
            //cout << "Assigned name " << name << endl; 
    }     
}

void Parser::load_children(MCNode* parent)
{
    MCNode* child;
    bool keep_reading = true;

    do
    {
        //cout << "adding child to " << parent->name << endl;
        child = parent->addChild();        
        load_node(child);

        switch (*character)
        {
            case ',':
                keep_reading = true;
                ++character;
                break;
            case ')':
                keep_reading = false;
                break;
        }
    }
    while (keep_reading);
    //cout << "terminating" << endl;
}

