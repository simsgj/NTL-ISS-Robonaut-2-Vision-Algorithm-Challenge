#ifndef GRAPH_H
#define GRAPH_H

#include "edge.h"
#include <algorithm>

using namespace std;

class Graph
{
public:
    Graph();
    void setEdges(Edge * edges, size_t nbEges);
    ~Graph();
    void sortEdges();
    Edge * getEdges();
    size_t getNbEdges();

protected:
    Edge * edges;
    size_t nbEdges;
};

#endif // GRAPH_H
