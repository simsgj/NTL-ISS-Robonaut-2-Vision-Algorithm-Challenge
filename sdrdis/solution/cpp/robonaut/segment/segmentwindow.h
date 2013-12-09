#ifndef SEGMENTWINDOW_H
#define SEGMENTWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include "image.h"
#include <QImage>
#include <QPixmap>
#include "imageutility.h"
#include "imagedisjointset.h"
#include "felzenhuttensegmentation.h"
#include <QRgb>
#include <imagesmoother.h>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QMap>
#include <QSlider>
#include <QMapIterator>
#include <QColor>
#include "segmentutility.h"
#include <QElapsedTimer>
#include "statistics.h"
#include "standardselector.h"
#include "ui_segmentwindow.h"
#include "felzenhuttensegmentation.h"
#include "standardmerger.h"

struct ValueProperties {
  float min;
  float max;
  int nbSteps;
  float value;
  QSlider * slider;
  QLabel * valueLabel;
};

namespace Ui {
class SegmentWindow;
}

class SegmentWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit SegmentWindow(QWidget *parent = 0);
  void setImage(Image<float> * image, bool (*isSetCorrect)(Set<float> *));
  ~SegmentWindow();
  void refresh();

  QPixmap getTerritoryMapPixmap(vector<int> * territoryMap);
  void addSliderToLayout(QVBoxLayout * layout, QString name, QString label, float min, float max, float value, int nbSteps = 100);
  QColor getGoldenRatioColor();
  void addMergeConfigurationToLayout(QVBoxLayout * layout, int index);
  void addSelectConfigurationToLayout(QVBoxLayout * layout, int index, float minSize = 0, float maxSize = 600, float minWidthHeight = 0.25, float maxWidthHeight = 4, float minSizeRatio = 0.15);

public slots:
  void setSliderValue(int newValue);
  
private:
  Ui::SegmentWindow *ui;
  Image<float> * image;
  QGraphicsScene * sceneFrom;
  QGraphicsScene * sceneSegmented;
  QGraphicsScene * sceneSelected1;
  QGraphicsScene * sceneMerged;
  QGraphicsScene * sceneSelected2;
  QGraphicsScene * sceneRectangles;
  QMap <QString, ValueProperties> settings;
  bool (*isSetCorrect)(Set<float> *);
};

#endif // SEGMENTWINDOW_H
