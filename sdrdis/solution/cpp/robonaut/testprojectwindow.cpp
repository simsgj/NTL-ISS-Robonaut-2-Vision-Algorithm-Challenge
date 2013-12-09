#include "testprojectwindow.h"
#include "ui_testprojectwindow.h"

TestProjectWindow::TestProjectWindow(QWidget *parent) :
  StandardImageViewerWindow(parent),
  ui(new Ui::TestProjectWindow)
{
  ui->setupUi(this);

  this->ui->image->installEventFilter(this);
  secondaryInfos = Utility::loadSecondaryObjects("data/model.csv", "data/secondary_positions/");
  this->selectImage(0);
}

TestProjectWindow::~TestProjectWindow()
{
  delete ui;
}


void TestProjectWindow::selectImage(int i) {
  if (i != selected) {
    positions.clear();
  }
  StandardImageViewerWindow::selectImage(i);

  vector<RealObjectInformations> imageSecondaryInfos = secondaryInfos.find(this->filesList.at(i / 2))->second;
  QPixmap leftImage(imagePath);
  vector<int> imageVector = Utility::loadImage(imagePath);
  Image<float> * image = RobonautEye::convertVectorToImage(imageVector);
  //QPixmap rightImage(folderPath + "RightImage_" + imageFilename);
  int fromWidth = leftImage.width();

  leftImage = leftImage.scaledToWidth(1024, Qt::SmoothTransformation);
  //rightImage = rightImage.scaledToWidth(512, Qt::SmoothTransformation);

  this->scale = leftImage.width() * 1.0 / fromWidth;

  QPainter p(&leftImage);
  p.setPen(QColor(255, 255, 255));

  vector <Point> primaryPositions = this->getPrimaryObjectsPositions();
  vector<RealObjectInformations> primaryObjectsRealInformations = this->getPrimaryObjectsRealInformations();
  primaryPositions.erase(primaryPositions.begin() + 1);
  primaryObjectsRealInformations.erase(primaryObjectsRealInformations.begin() + 1);
  Utility::drawPositionsWithIndex(p, primaryPositions, scale);

  p.setPen(QColor(0, 255, 0));
  Utility::drawPositions(p, positions, scale);



  vector <Point> primaryFilteredPositionsFrom;
  vector <Point> primaryFilteredPositionsTo;
  for (size_t j = 0; j < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; j++) {
    Point pointTo = primaryPositions.at(j);
    if (pointTo.x < 0 || pointTo.y < 0) {
      continue;
    }

    primaryFilteredPositionsFrom.push_back(TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position);
    primaryFilteredPositionsTo.push_back(primaryPositions.at(j));
  }

  if (primaryFilteredPositionsFrom.size() > 3) {

    Array2D<double> * H = Geometry::getPerspectiveProjectionHomography(
      &primaryFilteredPositionsFrom[0],
      &primaryFilteredPositionsTo[0],
      primaryFilteredPositionsFrom.size()
    );

    Point borders[4] = {
      Point(20, 19),
      Point(537, 19),
      Point(537, 770),
      Point(20, 770)
    };

    Point * projectedBorders = Geometry::getPerspectiveProjectionPoints(H, &borders[0], 4);

    p.setPen(QColor(255, 255, 255));
    Utility::drawPolygon(p, projectedBorders, 4, scale);

    delete projectedBorders;

    Point * projectionsRaw = Geometry::getPerspectiveProjectionPoints(H, RobonautEye().getSecondaryPoints(), 20);

    vector < Point > projections;
    for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
      projections.push_back(projectionsRaw[i]);
    }

    p.setPen(QColor(255, 0, 0));
    vector <Point> assignments = Utility::assignSecondaryInformations(imageSecondaryInfos, projections, i % 2 == 0);
    Utility::drawPositions(p, assignments, scale);



    p.setPen(QColor(0, 255, 255));

    Utility::drawPositions(p, projections, scale);
    RobonautEye eye;

    for (size_t j = 0; j < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; j++) {
      /*Point objectBorders[4] = {
        Point(
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.x - TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2,
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.y - TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2
        ),
        Point(
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.x + TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2,
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.y - TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2
        ),
        Point(
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.x + TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2,
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.y + TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2
        ),
        Point(
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.x - TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2,
          TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position.y + TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].size / 2
        )
      };

      Point * projectedObjectBorders = Geometry::getPerspectiveProjectionPoints(H, &objectBorders[0], 4);
      Utility::drawPolygon(p, projectedObjectBorders, 4, scale);
      delete projectedObjectBorders;
*/
/*primaryObjectsRealInformations.at(j).state == "HIDDEN"*/


      PrimaryObject primaryObjectProjection = RobonautEye::getPerspectiveProjection(H, TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j]);

      eye.computePrimaryObjectState(image, primaryObjectProjection);
      if (primaryObjectProjection.stateScores[0] > 0.5) {
        p.setPen(QColor(255, 0, 0));
      } else if (primaryObjectProjection.stateScores[1] > 0.5) {
        p.setPen(QColor(0, 255, 0));
      } else {
        p.setPen(QColor(255, 255, 0));
      }

      Utility::drawPrimaryObject(p, primaryObjectProjection, scale);
    }



  }
  delete image;

  /*
  if (positions.size() == 3) {
    OrthographicProjectionEquation eq = Geometry::getOrthographicProjectionEquation(
      TASKBOARD_SECONDARY_OBJECTS_POSITIONS[0],
      TASKBOARD_SECONDARY_OBJECTS_POSITIONS[4],
      TASKBOARD_SECONDARY_OBJECTS_POSITIONS[9],
      positions.at(0),
      positions.at(1),
      positions.at(2)
    );

    qDebug() << eq.xEq.a << eq.xEq.b << eq.xEq.c;
    qDebug() << eq.yEq.a << eq.yEq.b << eq.yEq.c;

    Point * projectionsRaw = Geometry::getOrthographicProjectionPoints(eq, &TASKBOARD_SECONDARY_OBJECTS_POSITIONS[0], TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS);

    vector < Point > projections;
    for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
      projections.push_back(projectionsRaw[i]);
    }

    p.setPen(QColor(255, 255, 0));

    Utility::drawPositions(p, projections, scale);
  }
  */

  if (positions.size() > 3) {
    size_t nb = positions.size();
    Point * pointsTo = new Point[nb];
    for (size_t i = 0; i < nb; i++) {
      pointsTo[i] = positions.at(i);
    }

    Array2D<double> * H = Geometry::getPerspectiveProjectionHomography(
      RobonautEye().getSecondaryPoints(),
      pointsTo,
      nb
    );

    Point * projectionsRaw = Geometry::getPerspectiveProjectionPoints(H, RobonautEye().getSecondaryPoints(), 20);

    vector < Point > projections;
    for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
      projections.push_back(projectionsRaw[i]);
    }

    p.setPen(QColor(255, 255, 0));

    Utility::drawPositions(p, projections, scale);
    delete pointsTo;
  }

  this->ui->image->setPixmap(leftImage);
  //this->ui->rightImage->setPixmap(rightImage);

}


bool TestProjectWindow::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::KeyRelease) {
    QKeyEvent * keyEvent = static_cast <QKeyEvent *> (event);
    if (keyEvent->key() == Qt::Key_Left) {
      this->selectImage(selected - 1);
    }
    if (keyEvent->key() == Qt::Key_Right) {
      this->selectImage(selected + 1);
    }
  }
  if (obj->objectName() == "image" && (event->type() == QEvent::MouseButtonRelease)) {
    QMouseEvent * mouseEvent = static_cast <QMouseEvent *> (event);

    if (mouseEvent->button() == 1) {
      this->addPoint(mouseEvent->x(), mouseEvent->y());
    } else {
      this->removePoint(mouseEvent->x(), mouseEvent->y());
    }
  }
  return false;
}

void TestProjectWindow::addPoint(int x, int y) {
  positions.push_back(Point(x / scale, y / scale));
  this->selectImage(selected);
}

void TestProjectWindow::removePoint(int x, int y) {
  float nearestDistance = -1;
  size_t nearest = 0;

  int realX = x / scale;
  int realY = y / scale;
  for (size_t i = 0; i < positions.size(); i++) {
    float distance = sqrt((realX - positions.at(i).x) * (realX - positions.at(i).x) + (realY - positions.at(i).y) * (realY - positions.at(i).y));
    if (nearestDistance < 0 || nearestDistance > distance) {
      nearest = i;
      nearestDistance = distance;
    }
  }

  if (nearestDistance > 0) {
    positions.erase(positions.begin() + nearest);
  }
  this->selectImage(selected);
}

