#ifndef FELZENHUTTENSEGMENTATION_H
#define FELZENHUTTENSEGMENTATION_H

#include "image.h"
#include "imagedisjointset.h"
#include "graph.h"
#include "edge.h"
#include <vector>

#define FELZEN_HUTTEN_THRESHOLD (size, c) (c/size)

/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/
// This code have been adapted in order to use this image class.

template <class T>
class FelzenHuttenSegmentation
{
public:
  FelzenHuttenSegmentation(float threshold, int minSize);
  Graph * getGraphFromImage(Image<T>* im);
  ImageDisjointSet<T> * segmentImage(Image<T> *im); // Image is considered as already smoothed
  ImageDisjointSet<T> * segmentGraph(int nbVertices, Graph * graph);
protected:
  float threshold;
  int minSize;
};

template <class T>
FelzenHuttenSegmentation<T>::FelzenHuttenSegmentation(float threshold, int minSize)
{
  this->threshold = threshold;
  this->minSize = minSize;
}

template <class T>
Graph * FelzenHuttenSegmentation<T>::getGraphFromImage(Image<T> *im) {
  int width = im->getWidth();
  int height = im->getHeight();

  // build graph
  size_t nbEdges = ((width-1)*(height-2)*4 + (height - 1) + (width-1) * 5);
  Edge * edges = new Edge[nbEdges];
  int num = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (x < width-1) {
        edges[num].from = y * width + x;
        edges[num].to = y * width + (x+1);
        edges[num].weight = im->diff(x, y, x+1, y);
        num++;
      }

      if (y < height-1) {
        edges[num].from = y * width + x;
        edges[num].to = (y+1) * width + x;
        edges[num].weight = im->diff(x, y, x, y+1);
        num++;
      }

      if ((x < width-1) && (y < height-1)) {
        edges[num].from = y * width + x;
        edges[num].to = (y+1) * width + (x+1);
        edges[num].weight = im->diff(x, y, x+1, y+1);
        num++;
      }

      if ((x < width-1) && (y > 0)) {
        edges[num].from = y * width + x;
        edges[num].to = (y-1) * width + (x+1);
        edges[num].weight = im->diff(x, y, x+1, y-1);
        num++;
      }
    }
  }

  Graph * graph = new Graph();
  graph->setEdges(edges, nbEdges);
  return graph;
}

template <class T>
ImageDisjointSet<T> * FelzenHuttenSegmentation<T>::segmentImage(Image<T> *im) {

  Graph * graph = this->getGraphFromImage(im);
  ImageDisjointSet<T> * disjointSet = this->segmentGraph(im->getWidth() * im->getHeight(), graph);

  Edge * edges = graph->getEdges();
  size_t nbEdges = graph->getNbEdges();
  // post process small components
  if (this->minSize > 0) {
    for (size_t i = 0; i < nbEdges; i++) {
      int a = disjointSet->find(edges[i].from);
      int b = disjointSet->find(edges[i].to);
      if ((a != b) && ((disjointSet->size(a) < this->minSize) || (disjointSet->size(b) < this->minSize))) {
        disjointSet->join(a, b);
      }
    }
  }

  delete graph;

  return disjointSet;


  //universe *u = segment_graph(width*height, num, edges, c);

  /*// segment


  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  delete [] edges;
  *num_ccs = u->num_sets();

  image<rgb> *output = new image<rgb>(width, height);

  // pick random colors for each component
  rgb *colors = new rgb[width*height];
  for (int i = 0; i < width*height; i++)
    colors[i] = random_rgb();

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int comp = u->find(y * width + x);
      imRef(output, x, y) = colors[comp];
    }
  }

  delete [] colors;
  delete u;
*/
}

template <class T>
ImageDisjointSet<T> * FelzenHuttenSegmentation<T>::segmentGraph(int nbVertices, Graph * graph) {
  graph->sortEdges();
  Edge * edges = graph->getEdges();
  size_t nbEdges = graph->getNbEdges();

  ImageDisjointSet<T> * disjointSet = new ImageDisjointSet<T>(nbVertices);


  // init thresholds
  float *thresholds = new float[nbVertices];
  for (int i = 0; i < nbVertices; i++) {
    thresholds[i] = this->threshold;
  }

  // for each edge, in non-decreasing weight order...
  for (size_t i = 0; i < nbEdges; i++) {
    Edge * currentEdge = &edges[i];

    // components conected by this edge
    int a = disjointSet->find(currentEdge->from);
    int b = disjointSet->find(currentEdge->to);
    if (a != b) {
      if ((currentEdge->weight <= thresholds[a]) && (currentEdge->weight <= thresholds[b])) {
        disjointSet->join(a, b);
        a = disjointSet->find(a);
        thresholds[a] = currentEdge->weight + (this->threshold / disjointSet->size(a));
      }
    }
  }

  // free up
  delete thresholds;


  return disjointSet;
  /*
  std::sort(edges, edges + num_edges);

  // make a disjoint-set forest
  universe *u = new universe(num_vertices);

  // init thresholds
  float *threshold = new float[num_vertices];
  for (int i = 0; i < num_vertices; i++)
    threshold[i] = THRESHOLD(1,c);

  // for each edge, in non-decreasing weight order...
  for (int i = 0; i < num_edges; i++) {
    edge *pedge = &edges[i];

    // components conected by this edge
    int a = u->find(pedge->a);
    int b = u->find(pedge->b);
    if (a != b) {
      if ((pedge->w <= threshold[a]) &&
    (pedge->w <= threshold[b])) {
  u->join(a, b);
  a = u->find(a);
  threshold[a] = pedge->w + THRESHOLD(u->size(a), c);
      }
    }
  }

  // free up
  delete threshold;
  return u;*/
}


#endif // FELZENHUTTENSEGMENTATION_H
