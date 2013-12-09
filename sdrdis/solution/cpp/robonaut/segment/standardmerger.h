#ifndef STANDARDMERGER_H
#define STANDARDMERGER_H


#include "image.h"
#include "imagedisjointset.h"
#include <vector>
#include <map>
#include <set>
#include "mergepair.h"

using namespace std;

/* Deleted but can be ideas
  min/max width
  min/max height
  min/max size ratio
  max gravity center distance
*/
struct StandardMergerSettings {
  float maxColorDiff;
  float minSharedThreshold;
  float coefSizeScore;
  float coefSharedThresholdScore;
  float coefColorDiffScore;
};



template <class T>
class StandardMerger
{
public:
  StandardMerger(StandardMergerSettings settings, int width, int height);
  void apply(vector< Set<T> * > * sets);
  bool canMergeSets(vector< Set<T> * > * sets, size_t first, size_t second);
  void mergeTwoSets(vector< Set<T> * > * sets, multimap<float, MergePair > * mergePairs, bool looseMerge);
  void mergeTwoSets(vector< Set<T> * > * sets, size_t first, size_t second);
  void addMergePair(multimap<float, MergePair > * mergePairs, vector< Set<T> * > * sets, size_t first, size_t second);
  void mergeAllSets(vector< Set<T> * > * sets, bool looseMerge);

protected:
  StandardMergerSettings settings;
  int width;
  int height;
};

template <class T>
StandardMerger<T>::StandardMerger(StandardMergerSettings settings, int width, int height)
{
  this->settings = settings;
  this->width = width;
  this->height = height;
}

template <class T>
void StandardMerger<T>::apply(vector< Set<T> * > * sets)
{
  this->mergeAllSets(sets, false);
}


template <class T>
void StandardMerger<T>::mergeAllSets(vector< Set<T> * > * sets, bool looseMerge) {
  multimap<float, MergePair > * mergePairs = new multimap<float, MergePair>();
  for (size_t i = 0; i < sets->size(); i++) {
    Set<T> * firstSet = sets->at(i);
    for (size_t j = i + 1; j < sets->size(); j++) {
      if (i != j) {
        if (firstSet->sharedThresholds[j] > 0 &&
            (looseMerge || this->canMergeSets(sets, i, j))
            ) {
          this->addMergePair(mergePairs, sets, i, j);
        }
      }
    }
  }


  while(mergePairs->size() > 0) {
    this->mergeTwoSets(sets, mergePairs, looseMerge);
  }

  delete mergePairs;
}

template <class T>
bool StandardMerger<T>::canMergeSets(vector< Set<T> * > * sets, size_t first, size_t second) {
  Set<T> * firstSet = sets->at(first);
  Set<T> * secondSet = sets->at(second);

  float sharedThreshold = fmax(
    firstSet->sharedThresholds[second] * 1.0 / firstSet->nbSharedThresholds,
    secondSet->sharedThresholds[first] * 1.0 / secondSet->nbSharedThresholds
  );


  if (fabs(firstSet->meanColor - secondSet->meanColor) > settings.maxColorDiff) {
    return false;
  }


  if (firstSet->sharedThresholds[second] == 0 || sharedThreshold < settings.minSharedThreshold) { // add
    return false;
  }

  return true;
}

template <class T>
void StandardMerger<T>::mergeTwoSets(vector< Set<T> * > * sets, multimap<float, MergePair > * mergePairs, bool looseMerge)
{
  MergePair mergePair = (*mergePairs->begin()).second;
  size_t first = mergePair.first;
  size_t second = mergePair.second;

  this->mergeTwoSets(sets, first, second);

  // Clear merge pairs with first or second Set... Second set doesn't exists anymore, and we will refresh merge pairs for first set
  for (multimap<float,  MergePair >::iterator it=mergePairs->begin(); it!=mergePairs->end();) {
    MergePair mergePair = it->second;
    if (mergePair.first == first || mergePair.second == first || mergePair.first == second || mergePair.second == second) {
      mergePairs->erase(it++);
      continue;
    }
    if (mergePair.first > second || mergePair.second > second) {
      if (mergePair.first > second) {
        mergePair.first--;
      }
      if (mergePair.second > second) {
        mergePair.second--;
      }
      it->second = mergePair;
    }

    ++it;
  }

  if (first > second) {
    first--;
  }
  Set<T> * firstSet = sets->at(first);

  for (size_t i = 0; i < sets->size(); i++) {
    if (i != first) {
      if (firstSet->sharedThresholds.at(i) > 0 &&
          (looseMerge || this->canMergeSets(sets, first, i))
          ) {
        this->addMergePair(mergePairs, sets, first, i);
      }
    }
  }
}

template <class T>
void StandardMerger<T>::mergeTwoSets(vector< Set<T> * > * sets, size_t first, size_t second)
{
  Set<T> * firstSet = sets->at(first);
  Set<T> * secondSet = sets->at(second);

  // Nb shared threshold isn't changed for any sets except firstSet since we are merging it with secondSet
  firstSet->nbSharedThresholds -= firstSet->sharedThresholds[second];

  // Merge shared threshold on other sets (into i)
  for (size_t i = 0; i < sets->size(); i++) {
    if (i != second) {
      Set<T> * currentSet = sets->at(i);
      int sharedThresholdWithSecond = currentSet->sharedThresholds[second];
      currentSet->sharedThresholds[first] += sharedThresholdWithSecond;
      if (i!= first) {
        firstSet->sharedThresholds[i] += sharedThresholdWithSecond;
        firstSet->nbSharedThresholds += sharedThresholdWithSecond;
      }
    }
  }

  // Erase shared threshold of second Set since it will not exist anymore
  for (size_t i = 0; i < sets->size(); i++) {
    Set<T> * currentSet = sets->at(i);
    currentSet->sharedThresholds.erase(currentSet->sharedThresholds.begin() + second);
  }

  // Merge two sets
  firstSet->minX = min(firstSet->minX, secondSet->minX);
  firstSet->maxX = max(firstSet->maxX, secondSet->maxX);
  firstSet->minY = min(firstSet->minY, secondSet->minY);
  firstSet->maxY = max(firstSet->maxY, secondSet->maxY);
  firstSet->gravityCenterX = (firstSet->gravityCenterX * firstSet->size + secondSet->gravityCenterX * secondSet->size) / (1.0 * firstSet->size + secondSet->size);
  firstSet->gravityCenterY = (firstSet->gravityCenterY * firstSet->size + secondSet->gravityCenterY * secondSet->size) / (1.0 * firstSet->size + secondSet->size);
  firstSet->size = firstSet->size + secondSet->size;
  firstSet->points.insert(firstSet->points.end(), secondSet->points.begin(), secondSet->points.end());
  firstSet->meanColor = (firstSet->meanColor * firstSet->size + secondSet->meanColor * secondSet->size) / (1.0 * firstSet->size + secondSet->size);

  // Removing second set
  sets->erase(sets->begin() + second);
  delete secondSet; // not sure
}

template <class T>
void StandardMerger<T>::addMergePair(multimap<float, MergePair > * mergePairs, vector< Set<T> * > * sets, size_t first, size_t second)
{
  Set<T> * firstSet = sets->at(first);
  Set<T> * secondSet = sets->at(second);

  float sizeScore = settings.coefSizeScore * fmax(firstSet->size, secondSet->size) / (width * height);
  float sharedThresholdScore = settings.coefSharedThresholdScore * (1.0 - fmax(
        firstSet->sharedThresholds[second] / firstSet->nbSharedThresholds,
        secondSet->sharedThresholds[first] / secondSet->nbSharedThresholds
      ));

  float colorDiffScore = settings.coefColorDiffScore * fabs(firstSet->meanColor - secondSet->meanColor) / 255.0;
  float sumScore = sqrt(
        sizeScore * sizeScore +
        sharedThresholdScore * sharedThresholdScore +
        colorDiffScore * colorDiffScore);

  MergePair mergePair;
  mergePair.first = first;
  mergePair.second = second;
  mergePair.score = sumScore;

  mergePairs->insert(pair<float, MergePair>(mergePair.score, mergePair));
}

#endif // STANDARDMERGER_H
