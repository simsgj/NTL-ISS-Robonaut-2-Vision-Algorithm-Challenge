#include "secondaryobjectclassifierwindow.h"
#include "ui_secondaryobjectclassifierwindow.h"

SecondaryObjectClassifierWindow::SecondaryObjectClassifierWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::SecondaryObjectClassifierWindow)
{
  ui->setupUi(this);
  this->installEventFilter(this);
  this->ui->image->installEventFilter(this);

  this->infos = Utility::loadCSV("data/model.csv");
  for (map<QString, vector <RealObjectInformations> >::iterator it=infos.begin(); it!=infos.end(); ++it) {
    filesList.push_back(it->first);
  }
}

SecondaryObjectClassifierWindow::~SecondaryObjectClassifierWindow()
{
  delete ui;
}


void SecondaryObjectClassifierWindow::selectImage(int i) {
  if (i < 0) {
    i = filesList.size() * 2 - 1;
  }
  if (i > (filesList.size() * 2 - 1)) { // after because unsigned
    i = 0;
  }

  selected = i;
  this->setWindowTitle(QString().setNum(i + 1) + " / " + QString().setNum(filesList.size() * 2));
  map<QString, vector <RealObjectInformations> >::iterator it = infos.find(filesList.at(i / 2)); //filesList.at(i / 2)
  QStringList pathList = it->first.split('/');
  pathList[1] = ((i % 2 == 0) ? "LeftImage_" : "RightImage_") + pathList[1];
  QString imagePath = "data/" + pathList.join('/');

  QImage im(imagePath);
  im = ImageUtility::changeGamma(im, 300);

  QPixmap rescaled = QPixmap::fromImage(im).scaledToWidth(1024, Qt::SmoothTransformation);

  scale = rescaled.width() * 1.0 / im.width();
  QPainter p(&rescaled);

  vector <int> imageVector = Utility::loadImage(imagePath);
  Image<float> * imageInstance = RobonautEye::convertVectorToImage(imageVector);
  RobonautEye eye;

  vector<ObjectInformations> squares;
  QElapsedTimer timer;
  timer.start();
  ObjectInformationsMapping * mapping = eye.getSecondaryObjectsInformationsFromImage(imageInstance);

  qDebug() << "ALIVE!" << mapping->positives.size();
  p.setPen(QColor(255, 0, 0));
  bool error;
  Array2D<double> * H = eye.getHomography(imageInstance, mapping, error);

  qDebug() << "ERROR:" << error;
  if (!error) {
    Point * primaryPoints = new Point[TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS];
    for (size_t j = 0; j < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; j++) {
      primaryPoints[j] = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position;
    }

    Point * projected = Geometry::getPerspectiveProjectionPoints(H, RobonautEye().getSecondaryPoints(), TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS);

    Utility::drawPositions(p, projected, TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS, scale);

    delete primaryPoints;

    delete projected;
  }

  delete H;
  delete imageInstance;

  statusBar()->showMessage("Processing time: " + QString().setNum(timer.elapsed()) + "ms");
  for (size_t i = 0; i < mapping->positives.size(); i++) {
    ObjectInformations positive = mapping->positives.at(i);
    positive.state = positive.secondaryObjectClassifierScore > SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD ? "1" : "0";
    squares.push_back(positive);
    p.drawText((positive.x * scale), ((positive.y + positive.height + 20) * scale), QString().setNum(i));
  }



  Utility::drawSquares(p, squares, scale);

  ui->image->setPixmap(rescaled);
}

bool SecondaryObjectClassifierWindow::eventFilter(QObject *obj, QEvent *event)
{
  if (obj->objectName() == "image" && (event->type() == QEvent::MouseButtonRelease)) {
    QMouseEvent * mouseEvent = static_cast <QMouseEvent *> (event);
    statusBar()->showMessage("Mouse position: " + QString().setNum(mouseEvent->x() / scale) + ", " + QString().setNum(mouseEvent->y() / scale));
  }
  if (event->type() == QEvent::KeyRelease) {
    QKeyEvent * keyEvent = static_cast <QKeyEvent *> (event);
    if (keyEvent->key() == Qt::Key_Left) {
      this->selectImage(selected - 1);
    }
    if (keyEvent->key() == Qt::Key_Right) {
      this->selectImage(selected + 1);
    }
  }
  return false;
}
