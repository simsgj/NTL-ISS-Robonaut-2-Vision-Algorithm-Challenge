#ifndef DISJOINTSET_H
#define DISJOINTSET_H

#include <map>
#include <vector>
#include "image.h"

using namespace std;

template <typename T>
struct Set {
  int minX;
  int maxX;
  int minY;
  int maxY;
  int gravityCenterX;
  int gravityCenterY;
  int size;
  vector <int> points;
  vector <int> sharedThresholds;
  int nbSharedThresholds;
  T meanColor;
};

typedef struct {
  int rank;
  int parent;
  // + firstAncestor ?
  int size;
} ds_elem;

template <class T>
class ImageDisjointSet
{
public:
  ImageDisjointSet(int nbElems);
  ~ImageDisjointSet();
  int find(int x);
  void join(int x, int y);
  int size(int x) const;
  int getNbSets() const;
  vector< Set<T> * > * getSets(Image<T> * im);
private:

protected:
  ds_elem *elems;
  int nbElems;
  int nbSets;
};

template <class T>
ImageDisjointSet<T>::ImageDisjointSet(int nbElems)
{
  this->elems = new ds_elem[nbElems];
  this->nbElems = nbElems;
  this->nbSets = nbElems;
  for (int i = 0; i < nbElems; i++) {
    elems[i].rank = 0;
    elems[i].size = 1;
    elems[i].parent = i;
  }
}

template <class T>
ImageDisjointSet<T>::~ImageDisjointSet() {
  delete [] elems;
}

template <class T>
int ImageDisjointSet<T>::find(int x) { // @todo might need some optimization
  int y = x;
  while (y != elems[y].parent) {
    y = elems[y].parent;
  }
  elems[x].parent = y;
  return y;
}

template <class T>
void ImageDisjointSet<T>::join(int x, int y) {
  if (elems[x].rank > elems[y].rank) {
    elems[y].parent = x;
    elems[x].size += elems[y].size;
  } else {
    elems[x].parent = y;
    elems[y].size += elems[x].size;
    if (elems[x].rank == elems[y].rank) {
      elems[y].rank++;
    }
  }
  nbSets--;
}

template <class T>
int ImageDisjointSet<T>::size(int x) const {
  return elems[x].size;
}

template <class T>
int ImageDisjointSet<T>::getNbSets() const {
  return nbSets;
}

template <class T>
vector< Set<T> * > * ImageDisjointSet<T>::getSets(Image<T> * im) {
  vector< Set<T> * > * sets = new vector< Set<T> * >();
  int width = im->getWidth();
  int height = im->getHeight();


  // drawing territory map
  vector<int> * territoryMap = new vector <int>(width * height);
  map<int, int> parentIdAssociation; // could use unordered_map here but since there is not C++ 11...
  map<int, int>::iterator it;
  int id = 0;
  int comp;
  int nbPoints;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      comp = this->find(y * width + x);
      it = parentIdAssociation.find(comp);
      if (it != parentIdAssociation.end()) {
        (*territoryMap)[x + y * width] = it->second;
        Set<T> * set = sets->at(it->second);
        if (x < set->minX) {
          set->minX = x;
        }
        if (x > set->maxX) {
          set->maxX = x;
        }
        if (y < set->minY) {
          set->minY = y;
        }
        if (y > set->maxY) {
          set->maxY = y;
        }
        nbPoints = set->points.size();
        set->gravityCenterX = (set->gravityCenterX * 1.0 * nbPoints + x) / (nbPoints + 1.0);
        set->gravityCenterY = (set->gravityCenterY * 1.0 * nbPoints + y) / (nbPoints + 1.0);
        set->meanColor = (set->meanColor * 1.0 * nbPoints + imgPixel(im, x, y)) / (nbPoints + 1.0);
        set->points.push_back(x + y * width);
        set->size++;
      } else {
        parentIdAssociation.insert(pair<int, int>(comp, id));
        (*territoryMap)[x + y * width] = id;
        Set<T> * set = new Set<T>;
        set->minX = x;
        set->maxX = x;
        set->minY = y;
        set->maxY = y;
        set->gravityCenterX = x;
        set->gravityCenterY = y;
        set->meanColor = imgPixel(im, x, y);
        set->points.push_back(x + y * width);
        set->nbSharedThresholds = 0;
        set->size = 1;
        sets->push_back(set);
        id++;
      }

    }
  }

  size_t nbSets = sets->size();
  for (size_t i = 0; i < nbSets; i++) {
    sets->at(i)->sharedThresholds.resize(nbSets, 0);
  }
  // Trier du plus petit au plus grand ou vice versa...

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int territoryFrom = (*territoryMap)[x + y * width];

      if (x > 0) {
        int leftTerritory = (*territoryMap)[x - 1 + y * width];
        sets->at(territoryFrom)->sharedThresholds[leftTerritory]++;
        if (leftTerritory != territoryFrom) {
          sets->at(territoryFrom)->nbSharedThresholds++;
        }
      }
      if (x < width - 1) {
        int rightTerritory = (*territoryMap)[x + 1 + y * width];
        sets->at(territoryFrom)->sharedThresholds[rightTerritory]++;
        if (rightTerritory != territoryFrom) {
          sets->at(territoryFrom)->nbSharedThresholds++;
        }
      }
      if (y > 0) {
        int topTerritory = (*territoryMap)[x + (y - 1) * width];
        sets->at(territoryFrom)->sharedThresholds[topTerritory]++;
        if (topTerritory != territoryFrom) {
          sets->at(territoryFrom)->nbSharedThresholds++;
        }
      }
      if (y < height - 1) {
        int bottomTerritory = (*territoryMap)[x + (y + 1) * width];
        sets->at(territoryFrom)->sharedThresholds[bottomTerritory]++;
        if (bottomTerritory != territoryFrom) {
          sets->at(territoryFrom)->nbSharedThresholds++;
        }
      }
    }
  }


  delete territoryMap;

  return sets;
}

#endif // DISJOINTSET_H
