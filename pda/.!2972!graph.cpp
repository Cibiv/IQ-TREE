/*
 * graph.cpp
 *
 *  Created on: Nov 14, 2013
 *      Author: olga
 */

#include "graph.h"
#include <iostream>
#include <list>
#include <limits.h>

Graph::Graph(int V){
    this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int v, int w){
