#include "segmentutility.h"

double SegmentUtility::goldenRatioH = 0;
double SegmentUtility::goldenRatioS = 0;

SegmentUtility::SegmentUtility()
{
}


QPixmap SegmentUtility::getTerritoryMapPixmap(vector<int> * territoryMap, int width, int height) {
  QImage im(width, height, QImage::Format_RGB32);

  SegmentUtility::goldenRatioH = 0;
  SegmentUtility::goldenRatioS = 0;
  QVector<QColor> colors;
  qDebug() << "SIZE BEFORE DISPLAY: " << territoryMap->size() << " " << width * height;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int comp = territoryMap->at(y * width + x);
      QColor color;

      if (comp > -1) {
        if (comp > colors.size() - 1) {
          int colorsSize = colors.size();
          for (int i = 0; i < (comp - colorsSize + 1); i++) {
            colors.push_back(getGoldenRatioColor());
          }
        }
        //qDebug() << "new color before " << comp << "/" << colors.size();
        color = colors[comp];
        //qDebug() << "new color after";
      } else {
        color.setRgb(255, 255, 255);
      }
      if (//(x > 0 && territoryMap->at(y * width + x) != territoryMap->at(y * width + x - 1)) ||
          (x < width - 1 && territoryMap->at(y * width + x) != territoryMap->at(y * width + x + 1)) ||
          //(y > 0 && territoryMap->at(y * width + x) != territoryMap->at((y - 1) * width + x)) ||
          (y < height - 1 && territoryMap->at(y * width + x) != territoryMap->at((y + 1) * width + x))
          ) {
        color = color.darker();
      }
      im.setPixel(x,y,color.rgb());
    }
  }


  return QPixmap::fromImage(im);
}

QColor SegmentUtility::getGoldenRatioColor() {
  SegmentUtility::goldenRatioH += 0.618033988749895;
  SegmentUtility::goldenRatioS += 0.74215;
  if (SegmentUtility::goldenRatioH > 1) {
    SegmentUtility::goldenRatioH -= 1;
  }
  if (SegmentUtility::goldenRatioS > 1) {
    SegmentUtility::goldenRatioS -= 1;
  }
  QColor color;
  color.setHsv(SegmentUtility::goldenRatioH * 255, 128 + (128 * SegmentUtility::goldenRatioS) - 64, 243 - 64 * SegmentUtility::goldenRatioS);
  return color;
}

QPixmap SegmentUtility::getSetsPixmap(vector< Set<float> * > * sets, int width, int height, bool showRectangles) {
  vector <int> * territoryMap = getTerritoryMapFromSets(sets, width, height);
  QPixmap pixmap = getTerritoryMapPixmap(territoryMap, width, height);
  if (showRectangles) {
    QPainter p(&pixmap);
    p.setPen(QColor().black());
    for (size_t i = 0; i < sets->size(); i++) {
      Set <float> * currentSet = sets->at(i);
      p.drawRect(currentSet->minX, currentSet->minY, currentSet->maxX - currentSet->minX, currentSet->maxY - currentSet->minY);
    }
  }
  delete territoryMap;
  return pixmap;
}

vector <int> * SegmentUtility::getTerritoryMapFromSets(vector< Set<float> * > * sets, int width, int height) {
  vector <int> * territoryMap;
  territoryMap = new vector <int>();
  territoryMap->resize(width * height, -1);
  for (size_t i = 0; i < sets->size(); i++) {
    vector <int> points = sets->at(i)->points;
    for (vector<int>::iterator it=points.begin(); it < points.end(); ++it) {
      (*territoryMap)[(*it)] = i;
    }
  }
  return territoryMap;
}
