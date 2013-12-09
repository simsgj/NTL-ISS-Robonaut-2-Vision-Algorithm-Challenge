#include "segmentwindow.h"


SegmentWindow::SegmentWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::SegmentWindow)
{
  ui->setupUi(this);

  this->sceneFrom = new QGraphicsScene(this);
  ui->viewFrom->setScene(this->sceneFrom);

  this->sceneSegmented = new QGraphicsScene(this);
  ui->viewSegmented->setScene(this->sceneSegmented);

  this->sceneSelected1 = new QGraphicsScene(this);
  ui->viewSelected1->setScene(this->sceneSelected1);

  this->sceneMerged = new QGraphicsScene(this);
  ui->viewMerged->setScene(this->sceneMerged);

  this->sceneSelected2 = new QGraphicsScene(this);
  ui->viewSelected2->setScene(this->sceneSelected2);

  this->sceneRectangles = new QGraphicsScene(this);
  ui->viewRectangles->setScene(this->sceneRectangles);

  ui->preprocessingLayout->setAlignment(Qt::AlignTop);
  addSliderToLayout(ui->preprocessingLayout, "blur", "Blur: ", 0.04, 4, 1.1, 198);
  addSliderToLayout(ui->preprocessingLayout, "sigma", "Contrast: ", 0, 128, 100, 128);

  ui->splitLayout->setAlignment(Qt::AlignTop);
  addSliderToLayout(ui->splitLayout, "splitThreshold", "Threshold: ", 10, 1000, 200, 990);
  addSliderToLayout(ui->splitLayout, "splitMinSize", "Min size: ", 0, 1000, 30, 1000);
/*
  firstSettings.minSize = 50;
  firstSettings.maxSize = 10000;
  firstSettings.minWidthHeightRatio = 0.5;
  firstSettings.maxWidthHeightRatio = 5;
  firstSettings.minSizeRatio = 0;
  firstSettings.maxSizeRatio = 1000;
  firstSettings.maxColorDiff = 128;
  firstSettings.minSharedThreshold = 0.1;
  firstSettings.maxCenterOfGravityDistance = 1000;
  firstSettings.coefSizeScore = 1;
  firstSettings.coefSharedThresholdScore = 1;
  firstSettings.coefColorDiffScore = 1;
  firstSettings.coefGravityCenterDistanceScore = 1;
*/
  ui->selectLayout1->setAlignment(Qt::AlignTop);
  this->addSelectConfigurationToLayout(ui->selectLayout1, 0);

  ui->mergeLayout->setAlignment(Qt::AlignTop);
  this->addMergeConfigurationToLayout(ui->mergeLayout, 0);

  ui->selectLayout2->setAlignment(Qt::AlignTop);
  this->addSelectConfigurationToLayout(ui->selectLayout2, 1, 0, 250, 0.5, 2, 0.5);

  /*ui->mergeLayout2->setAlignment(Qt::AlignTop);
  this->addMergeConfigurationToLayout(ui->mergeLayout2, 1);*/
}

void SegmentWindow::addMergeConfigurationToLayout(QVBoxLayout * layout, int index) {
  addSliderToLayout(layout, "maxColorDiff" + QString().setNum(index), "Max color diff.: ", 0, 255, 50, 255);
  addSliderToLayout(layout, "minSharedThreshold" + QString().setNum(index), "Min shared thresh.: ", 0, 1, 0.2, 100);
}

void SegmentWindow::addSelectConfigurationToLayout(QVBoxLayout * layout, int index, float minSize, float maxSize, float minWidthHeight, float maxWidthHeight, float minSizeRatio) {
  addSliderToLayout(layout, "minSize" + QString().setNum(index), "Min size: ", 0, 200, minSize, 200);
  addSliderToLayout(layout, "maxSize" + QString().setNum(index), "Max size: ", 0, 1000, maxSize, 1000);
  addSliderToLayout(layout, "minWidthHeightRatio" + QString().setNum(index), "Min width/height: ", 0, 10, minWidthHeight, 1000);
  addSliderToLayout(layout, "maxWidthHeightRatio" + QString().setNum(index), "Max width/height: ", 0, 10, maxWidthHeight, 1000);
  addSliderToLayout(layout, "minSizeBoxRatio" + QString().setNum(index), "Min size/box: ", 0, 1, minSizeRatio, 100);
}

void SegmentWindow::addSliderToLayout(QVBoxLayout * layout, QString id, QString label, float min, float max, float value, int nbSteps) {
  QHBoxLayout * horizontalLayout = new QHBoxLayout();


  QLabel * nameLabel = new QLabel(label);
  nameLabel->setMinimumWidth(50);

  QSlider * slider = new QSlider(Qt::Horizontal);
  slider->setMaximum(nbSteps);
  slider->setValue(((value - min) / (max - min)) * nbSteps);

  QLabel * valueLabel = new QLabel(QString().setNum(value));
  valueLabel->setMinimumWidth(50);

  horizontalLayout->addWidget(nameLabel);
  horizontalLayout->addWidget(slider);
  horizontalLayout->addWidget(valueLabel);
  layout->addLayout(horizontalLayout);

  ValueProperties prop;
  prop.min = min;
  prop.max = max;
  prop.value = value;
  prop.nbSteps = nbSteps;
  prop.slider = slider;
  prop.valueLabel = valueLabel;

  settings.insert(id, prop);

  QObject::connect(slider, SIGNAL(valueChanged(int)),
                        this, SLOT(setSliderValue(int)));
}

void SegmentWindow::setSliderValue(int newValue) {
  QMapIterator<QString, ValueProperties > i(settings);
  while (i.hasNext()) {
      i.next();
      ValueProperties prop = i.value();
      prop.value = (prop.slider->value() * 1.0 / prop.nbSteps) * (prop.max - prop.min) + prop.min;
      prop.valueLabel->setText(QString().setNum(prop.value));
      settings[i.key()] = prop;
  }
  this->refresh();
}

SegmentWindow::~SegmentWindow()
{
  delete ui;
}

void SegmentWindow::setImage(Image<float> * image, bool (*isSetCorrect)(Set<float> *)) {
    this->image = image;
    this->isSetCorrect = isSetCorrect;
    this->refresh();
}

void SegmentWindow::refresh() {
  statusBar()->showMessage("In progress...");
  QElapsedTimer overallTimer;
  overallTimer.start();

  QElapsedTimer timer;
  qint32 smoothTime;
  qint32 segmentTime;
  qint32 disjointToSetTime;
  qint32 mergeTime;
  qint32 select1Time;
  qint32 select2Time;

  timer.start();
  Image<float> * im = ImageSmoother<float>::apply(this->image, settings["blur"].value);
  Statistics<float>::normalize(im->img, im->getWidth() * im->getHeight(), 128, settings["sigma"].value);
  Statistics<float>::applyLowerBound(im->img, im->getWidth() * im->getHeight(), 0);
  Statistics<float>::applyUpperBound(im->img, im->getWidth() * im->getHeight(), 255);
  smoothTime = timer.elapsed();
  QPixmap pixmapFrom = ImageUtility::imageFloatToQPixmap(im);

  timer.restart();
  FelzenHuttenSegmentation<float> segmenter(settings["splitThreshold"].value, settings["splitMinSize"].value);
  ImageDisjointSet<float> * disjointSet = segmenter.segmentImage(im);
  segmentTime = timer.elapsed();
  //QPixmap pixmapSegmented = SegmentWindow::getDisjointSetPixmap(disjointSet);


  StandardSelectorSettings firstSelectSettings;
  firstSelectSettings.minSize = settings["minSize0"].value;
  firstSelectSettings.maxSize = settings["maxSize0"].value;
  firstSelectSettings.minWidthHeightRatio = settings["minWidthHeightRatio0"].value;
  firstSelectSettings.maxWidthHeightRatio = settings["maxWidthHeightRatio0"].value;
  firstSelectSettings.minSizeBoxRatio = settings["minSizeBoxRatio0"].value;
  StandardSelector<float> firstSelector(firstSelectSettings, im->getWidth(), im->getHeight());

  StandardSelectorSettings secondSelectSettings;
  secondSelectSettings.minSize = settings["minSize1"].value;
  secondSelectSettings.maxSize = settings["maxSize1"].value;
  secondSelectSettings.minWidthHeightRatio = settings["minWidthHeightRatio1"].value;
  secondSelectSettings.maxWidthHeightRatio = settings["maxWidthHeightRatio1"].value;
  secondSelectSettings.minSizeBoxRatio = settings["minSizeBoxRatio1"].value;
  StandardSelector<float> secondSelector(secondSelectSettings, im->getWidth(), im->getHeight());


  StandardMergerSettings mergeSettings;
  mergeSettings.maxColorDiff = settings["maxColorDiff0"].value;
  mergeSettings.minSharedThreshold = settings["minSharedThreshold0"].value;
  mergeSettings.coefSizeScore = 1;
  mergeSettings.coefSharedThresholdScore = 1;
  mergeSettings.coefColorDiffScore = 1;
  StandardMerger<float> merger(mergeSettings, im->getWidth(), im->getHeight());

  timer.restart();
  vector< Set<float> * > * sets = disjointSet->getSets(im);
  disjointToSetTime = timer.elapsed();
  QPixmap pixmapSegmented = SegmentUtility::getSetsPixmap(sets, im->getWidth(), im->getHeight());

  timer.restart();
  firstSelector.apply(sets);
  select1Time = timer.elapsed();
  QPixmap pixmapSelected1 = SegmentUtility::getSetsPixmap(sets, im->getWidth(), im->getHeight());

  timer.restart();
  merger.apply(sets);
  mergeTime = timer.elapsed();
  QPixmap pixmapMerged = SegmentUtility::getSetsPixmap(sets, im->getWidth(), im->getHeight());

  timer.restart();
  secondSelector.apply(sets);
  select2Time = timer.elapsed();
  QPixmap pixmapSelected2 = SegmentUtility::getSetsPixmap(sets, im->getWidth(), im->getHeight(), true);

  QPixmap pixmapRectangles = pixmapFrom;
  QPainter p(&pixmapRectangles);

  for (size_t i = 0; i < sets->size(); i++) {
    Set<float> * currentSet = sets->at(i);
    if (this->isSetCorrect(currentSet)) {
      p.setPen(QColor(0, 255, 0));
    } else {
      p.setPen(QColor(255, 0, 0));
    }
    p.drawRect(currentSet->minX, currentSet->minY, currentSet->maxX - currentSet->minX, currentSet->maxY - currentSet->minY);
  }

  delete disjointSet;
  delete sets;
  delete im;

  this->sceneFrom->clear();
  this->sceneFrom->addPixmap(pixmapFrom);
  this->sceneFrom->setSceneRect(0, 0, pixmapFrom.width(), pixmapFrom.height());

  /*
  this->sceneSegmented->clear();
  this->sceneSegmented->addPixmap(pixmapSegmented);
  this->sceneSegmented->setSceneRect(0, 0, pixmapSegmented.width(), pixmapSegmented.height());
*/

  this->sceneSegmented->clear();
  this->sceneSegmented->addPixmap(pixmapSegmented);
  this->sceneSegmented->setSceneRect(0, 0, pixmapSegmented.width(), pixmapSegmented.height());

  this->sceneSelected1->clear();
  this->sceneSelected1->addPixmap(pixmapSelected1);
  this->sceneSelected1->setSceneRect(0, 0, pixmapSelected1.width(), pixmapSelected1.height());

  this->sceneMerged->clear();
  this->sceneMerged->addPixmap(pixmapMerged);
  this->sceneMerged->setSceneRect(0, 0, pixmapMerged.width(), pixmapMerged.height());

  this->sceneSelected2->clear();
  this->sceneSelected2->addPixmap(pixmapSelected2);
  this->sceneSelected2->setSceneRect(0, 0, pixmapSelected2.width(), pixmapSelected2.height());

  this->sceneRectangles->clear();
  this->sceneRectangles->addPixmap(pixmapRectangles);
  this->sceneRectangles->setSceneRect(0, 0, pixmapRectangles.width(), pixmapRectangles.height());

  qint64 overall = overallTimer.elapsed();
  statusBar()->showMessage(
        QString("Overall process time: ") + QString().setNum(overall) + QString("ms | ") +
        QString("Smooth: ") + QString().setNum(smoothTime) + QString("ms | ") +
        QString("Segment: ") + QString().setNum(segmentTime) + QString("ms | ") +
        QString("Select/Merge: ") + QString().setNum(disjointToSetTime + mergeTime + select1Time) + QString("ms (") + QString().setNum(disjointToSetTime) + QString(", ") + QString().setNum(select1Time) + QString(", ") + QString().setNum(mergeTime) + QString(", ") + QString().setNum(select2Time) + QString(") | ") +
        QString("Remaining: ") + QString().setNum(overall - smoothTime - segmentTime - disjointToSetTime - select1Time - mergeTime - select2Time) + QString("ms")
        );
}


