#include "overallwindow.h"
#include "ui_overallwindow.h"

OverallWindow::OverallWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::OverallWindow)
{
  ui->setupUi(this);
  this->installEventFilter(this);
  loadInformations();
}

OverallWindow::~OverallWindow()
{
  delete ui;
}

void OverallWindow::loadInformations() {
  this->infos = Utility::loadCSV("data/model.csv");
  for (map<QString, vector <RealObjectInformations> >::iterator it=infos.begin(); it!=infos.end(); ++it) {
    filesList.push_back(it->first);
  }
}

void OverallWindow::selectImage(int i) {
  if (i < 0) {
    i = filesList.size() - 1;
  }
  if (i > (filesList.size() - 1)) { // after because unsigned
    i = 0;
  }
  selected = i;
  this->setWindowTitle(QString().setNum(i + 1) + " / " + QString().setNum(filesList.size()));

  QStringList pathList = filesList.at(i).split('/');
  QString pathSecondItem = pathList[1];
  pathList[1] = "LeftImage_"+ pathSecondItem;
  QString leftImagePath = "data/" + pathList.join('/');
  pathList[1] = "RightImage_"+ pathSecondItem;
  QString rightImagePath = "data/" + pathList.join('/');

  QPixmap leftImage(leftImagePath);
  QPixmap rightImage(rightImagePath);

  QPixmap scaledLeftImage = leftImage.scaledToWidth(800, Qt::SmoothTransformation);
  QPixmap scaledRightImage = rightImage.scaledToWidth(800, Qt::SmoothTransformation);
  QPainter leftPainter(&scaledLeftImage);
  QPainter rightPainter(&scaledRightImage);

  float scale = scaledLeftImage.width()  * 1.0 / leftImage.width();
  /*drawPrimaryObjects(leftPainter, leftImagePath, scale);
  drawPrimaryObjects(rightPainter, rightImagePath, scale);*/
  drawFinalObjects(leftPainter, rightPainter, leftImagePath, rightImagePath, scale);


  this->ui->leftImage->setPixmap(scaledLeftImage);
  this->ui->rightImage->setPixmap(scaledRightImage);
}

void OverallWindow::drawFinalObjects(QPainter &leftPainter, QPainter &rightPainter, QString leftImagePath, QString rightImagePath, float scale) {
  vector<RealObjectInformations> infos = this->infos.find(filesList.at(selected))->second;
  vector<int> leftImageVector = Utility::loadImage(leftImagePath);
  Image<float> * leftImage = RobonautEye::convertVectorToImage(leftImageVector);
  vector<int> rightImageVector = Utility::loadImage(rightImagePath);
  Image<float> * rightImage = RobonautEye::convertVectorToImage(rightImageVector);
  RobonautEye eye;
  FinalObject * finalObjects = eye.getFinalObjects(leftImage, rightImage);
  float positionsScore = 0;
  float hiddenScore = 0;
  float stateScore = 0;
  float threshold = 10; // MUST BE CHANGED !
  qDebug() << "SWITCH FINAL POSITION:" << finalObjects[3].leftPosition.x << finalObjects[3].leftPosition.y << finalObjects[3].rightPosition.x << finalObjects[3].rightPosition.y;
  qDebug() << "POWER SWITCH FINAL POSITION:" << finalObjects[1].leftPosition.x << finalObjects[1].leftPosition.y << finalObjects[1].rightPosition.x << finalObjects[1].rightPosition.y;
  for (size_t i = 0; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    RealObjectInformations infosItem = infos.at(i);
    QColor color;
    if (finalObjects[i].state == 0) {
      color = QColor(255, 0, 0);
    } else if (finalObjects[i].state == 1) {
      color = QColor(0, 255, 0);
    } else {
      color = QColor(255, 255, 0);
    }
    leftPainter.setPen(color);
    rightPainter.setPen(color);
    int leftX = finalObjects[i].leftPosition.x * scale;
    int leftY = finalObjects[i].leftPosition.y * scale;
    int rightX = finalObjects[i].rightPosition.x * scale;
    int rightY = finalObjects[i].rightPosition.y * scale;
    leftPainter.drawLine(leftX - 3, leftY, leftX + 3, leftY);
    leftPainter.drawLine(leftX, leftY - 3, leftX, leftY + 3);
    rightPainter.drawLine(rightX - 3, rightY, rightX + 3, rightY);
    rightPainter.drawLine(rightX, rightY - 3, rightX, rightY + 3);

    int originalLeftX = infosItem.leftX * scale;
    int originalLeftY = infosItem.leftY * scale;
    int originalRightX = infosItem.rightX * scale;
    int originalRightY = infosItem.rightY * scale;

    leftPainter.setPen(QColor(255, 255, 255));
    rightPainter.setPen(QColor(255, 255, 255));

    leftPainter.drawLine(originalLeftX - 3, originalLeftY, originalLeftX + 3, originalLeftY);
    leftPainter.drawLine(originalLeftX, originalLeftY - 3, originalLeftX, originalLeftY + 3);
    rightPainter.drawLine(originalRightX - 3, originalRightY, originalRightX + 3, originalRightY);
    rightPainter.drawLine(originalRightX, originalRightY - 3, originalRightX, originalRightY + 3);

    if (infosItem.state == "HIDDEN") {
      hiddenScore += 3;
    } else {
      if (infosItem.leftX != -1 || infosItem.leftY != -1) {
        float distance = sqrt((finalObjects[i].leftPosition.x - infosItem.leftX) * (finalObjects[i].leftPosition.x - infosItem.leftX) + (finalObjects[i].leftPosition.y - infosItem.leftY) * (finalObjects[i].leftPosition.y - infosItem.leftY));
        if (distance <= 2 * threshold) {
          if (distance <= threshold) {
            positionsScore += 1;
          } else {
            positionsScore += 1 - pow(((distance - threshold) * 1.0 / threshold), 1.5);
          }
        }
      } else {
        positionsScore += 1;
      }

      if (infosItem.rightX != -1 || infosItem.rightY != -1) {
        float distance = sqrt((finalObjects[i].rightPosition.x - infosItem.rightX) * (finalObjects[i].rightPosition.x - infosItem.rightX) + (finalObjects[i].rightPosition.y - infosItem.rightY) * (finalObjects[i].rightPosition.y - infosItem.rightY));
        if (distance <= 2 * threshold) {
          if (distance <= threshold) {
            positionsScore += 1;
          } else {
            positionsScore += 1 - pow(((distance - threshold) * 1.0 / threshold), 1.5);
          }
        }
      } else {
        positionsScore += 1;
      }

      if (((infosItem.state == "DOWN" || infosItem.state == "OFF") && finalObjects[i].state == 0) ||
          ((infosItem.state == "UP" || infosItem.state == "ON") && finalObjects[i].state == 1) ||
          ((infosItem.state == "CENTER") && finalObjects[i].state == 2)
          ) {
        stateScore += 1;
      }
    }


  }

  qDebug() << "POSITION SCORE:" << positionsScore << "(" << ((positionsScore) * 1000000.0 / 66.0) << ")";
  qDebug() << "HIDDEN SCORE:" << hiddenScore << "(" << ((hiddenScore) * 1000000.0 / 66.0) << ")";
  qDebug() << "STATE SCORE:" << stateScore << "(" << ((stateScore) * 1000000.0 / 66.0) << ")";
  qDebug() << "TOTAL:" << positionsScore + hiddenScore + stateScore << "(" << ((positionsScore + hiddenScore + stateScore) * 1000000.0 / 66.0) << ")";

}

void OverallWindow::drawPrimaryObjects(QPainter &p, QString imagePath, float scale) {
  RobonautEye eye;

  vector<int> imageVector = Utility::loadImage(imagePath);
  Image<float> * imageInstance = RobonautEye::convertVectorToImage(imageVector);
  PrimaryObject * primaryObjects = eye.processImage(imageInstance);

  for (size_t i = 0; i < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; i++) {
    if (primaryObjects[i].stateScores[0] > 0.5) {
      p.setPen(QColor(255, 0, 0));
    } else if (primaryObjects[i].stateScores[1] > 0.5) {
      p.setPen(QColor(0, 255, 0));
    } else {
      p.setPen(QColor(255, 255, 0));
    }

    Utility::drawPrimaryObject(p, primaryObjects[i], scale);
  }

  delete imageInstance;
  delete primaryObjects;
}

bool OverallWindow::eventFilter(QObject *obj, QEvent *event)
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
  return false;
}
