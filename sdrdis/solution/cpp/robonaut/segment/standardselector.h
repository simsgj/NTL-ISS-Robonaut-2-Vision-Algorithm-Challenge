#ifndef STANDARDSELECTOR_H
#define STANDARDSELECTOR_H

#include "imagedisjointset.h"

struct StandardSelectorSettings {
  float minSize;
  float maxSize;
  float minWidthHeightRatio;
  float maxWidthHeightRatio;
  float minSizeBoxRatio;
};

template <class T>
class StandardSelector
{
public:
  StandardSelector(StandardSelectorSettings settings, int width, int height);
  bool isSetRespectingSettings(Set<T> * set, bool testRatio = false, bool testMinimums = false);
  void apply(vector< Set<T> * > * sets);
  void deleteSet(vector< Set<T> * > * sets, size_t index);

protected:
  StandardSelectorSettings settings;
  int width;
  int height;
};

template <class T>
StandardSelector<T>::StandardSelector(StandardSelectorSettings settings, int width, int height)
{
  this->settings = settings;
  this->width = width;
  this->height = height;
}


template <class T>
bool StandardSelector<T>::isSetRespectingSettings(Set<T> * set, bool testRatio, bool testMinimums) {
  if (set->size > settings.maxSize) {
    return false;
  }

  if (testRatio) {
    if ((set->maxY - set->minY) == 0) {
      return false;
    }
    float widthHeightRatio = (set->maxX - set->minX) * 1.0 / (set->maxY - set->minY);
    if (testRatio && (widthHeightRatio < settings.minWidthHeightRatio || widthHeightRatio > settings.maxWidthHeightRatio)) {
      return false;
    }

    float sizeBoxRatio = set->size * 1.0 / ((set->maxX - set->minX + 1) * (set->maxY - set->minY + 1));
    if (sizeBoxRatio < settings.minSizeBoxRatio) {
      return false;
    }
  }

  if (testMinimums) {
    if (set->size < settings.minSize) {
      return false;
    }
  }

  return true;
}


template <class T>
void StandardSelector<T>::apply(vector< Set<T> * > * sets) {
  for (size_t i = 0; i < sets->size(); i++) {
    Set<T> * currentSet = sets->at(i);
    if (!this->isSetRespectingSettings(currentSet, true, true)) {
      this->deleteSet(sets, i);
      i--;
    }
  }
}

template <class T>
void StandardSelector<T>::deleteSet(vector< Set<T> * > * sets, size_t index) {
  delete sets->at(index);
  sets->erase(sets->begin() + index);
  for (size_t i = 0; i < sets->size(); i++) {
    Set<T> * currentSet = sets->at(i);
    //currentSet->nbSharedThresholds -= currentSet->sharedThresholds[index]; -> commented because when we delete a set, we want the sharedThreshold ratio to stay the same on the others
    currentSet->sharedThresholds.erase(currentSet->sharedThresholds.begin() + index);
  }
}



#endif // STANDARDSELECTOR_H
