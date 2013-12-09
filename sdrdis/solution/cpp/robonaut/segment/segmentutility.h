#ifndef SEGMENTUTILITY_H
#define SEGMENTUTILITY_H

#include <QDebug>
#include <QColor>
#include <QImage>
#include <QPixmap>
#include <QPainter>
#include <vector>
#include "imagedisjointset.h"

using namespace std;

class SegmentUtility
{
public:
  SegmentUtility();
  static QPixmap getSetsPixmap(vector< Set<float> * > * sets, int width, int height, bool showRectangles = false);
  static vector <int> * getTerritoryMapFromSets(vector< Set<float> * > * sets, int width, int height);
  static QPixmap getTerritoryMapPixmap(vector<int> * territoryMap, int width, int height);
  static QColor getGoldenRatioColor();

protected:
  static double goldenRatioH;
  static double goldenRatioS;
};

#endif // SEGMENTUTILITY_H
