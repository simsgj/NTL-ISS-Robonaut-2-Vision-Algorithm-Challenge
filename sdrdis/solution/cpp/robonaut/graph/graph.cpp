#include "graph.h"

Graph::Graph()
{
}

Graph::~Graph() {
    delete this->edges;
}

void Graph::setEdges(Edge * edges, size_t nbEdges) {
  this->edges = edges;
  this->nbEdges = nbEdges;
}

Edge * Graph::getEdges() {
  return this->edges;
}

void Graph::sortEdges() {
  sort(this->edges, this->edges + this->nbEdges);
}

size_t Graph::getNbEdges() {
  return this->nbEdges;
}
