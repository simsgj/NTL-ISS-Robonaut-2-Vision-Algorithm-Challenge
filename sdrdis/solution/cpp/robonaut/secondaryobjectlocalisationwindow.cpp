#include "secondaryobjectlocalisationwindow.h"
#include "ui_secondaryobjectlocalisationwindow.h"

SecondaryObjectLocalisationWindow::SecondaryObjectLocalisationWindow(QWidget *parent, bool localisationMode) :
  QMainWindow(parent),
  ui(new Ui::SecondaryObjectLocalisationWindow)
{
  ui->setupUi(this);
  ui->image->installEventFilter(this);
  this->installEventFilter(this);
  this->localisationMode = localisationMode;
  if (this->localisationMode) {
    this->infos = Utility::loadCSV("data/model.csv");
  } else {
    this->infos = Utility::loadSecondaryObjects("data/model.csv", "data/secondary_positions/");
  }
  for (map<QString, vector <RealObjectInformations> >::iterator it=infos.begin(); it!=infos.end(); ++it) {
    filesList.push_back(it->first);
  }
  //selectImage(179);
}

SecondaryObjectLocalisationWindow::~SecondaryObjectLocalisationWindow()
{
  delete ui;
}

void SecondaryObjectLocalisationWindow::selectImage(int i) {
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
  imagePath = "data/" + pathList.join('/');

  QImage im(imagePath);
  im = ImageUtility::changeGamma(im, 300);

  QPixmap rescaled = QPixmap::fromImage(im).scaledToWidth(1024, Qt::SmoothTransformation);



  QPainter p(&rescaled);
  p.setPen(QColor(0, 255, 0));
  scale = rescaled.width() * 1.0 / im.width();
  Utility::drawInformations(p, it->second, i % 2 == 0, scale);

  if (this->localisationMode) {
    secondaryObjectsPath = "data/secondary_positions/" + pathList.join('_') + ".txt";
    loadSecondaryObjectsFile();
    p.setPen(QColor(255, 0, 0));
    Utility::drawSecondaryObjects(p, secondaryPositions, scale);
  } else {
    squarePath = "data/squares/" + pathList.join('_') + ".txt";
    if (QFile(squarePath).exists()) {
      this->squares = Utility::loadObjectsInformations(squarePath);

      /*
      // TEMPORARY IN ORDER TO MIGRATE TO SQUARES WITH LETTERS AND NUMBERS
      for (size_t j = 0; j < this->squares.size(); j++) {
        ObjectInformations currentObject = this->squares.at(j);
        int centerX = currentObject.x + currentObject.width / 2;
        int centerY = currentObject.y + currentObject.height / 2;
        for (size_t k = 0; k < it->second.size(); k++) {
          RealObjectInformations pos = it->second.at(k);
          if (pos.state == "TERTIARY") {
            int posX = i % 2 == 0 ? pos.leftX : pos.rightX;
            int posY = i % 2 == 0 ? pos.leftY : pos.rightY;
            if (sqrt((centerX - posX) * (centerX - posX) + (centerY - posY) * (centerY - posY)) < 15) {
              this->squares.at(j).state = "1";
              break;
            }
          }
        }
      }

      Utility::saveObjectsInformations(squarePath, this->squares);
      */
    } else {
      vector <int> imageVector = Utility::loadImage(imagePath);
      Image<float> * image = RobonautEye::convertVectorToImage(imageVector);
      RobonautEye eye;

      this->squares.clear();

      ObjectInformationsMapping mapping = *eye.getSecondaryObjectsInformationsFromImage(image);

      this->squares = mapping.positives;
/*
      for (map<int, ObjectInformations>::iterator it=mapping.objectInformationsMap.begin(); it!=mapping.objectInformationsMap.end(); ++it) {
        this->squares.push_back(it->second);
      }
*/
      delete image;


      for (size_t j = 0; j < this->squares.size(); j++) {
        this->squares.at(j).state = "0";
        ObjectInformations currentObject = this->squares.at(j);
        int centerX = currentObject.x + currentObject.width / 2;
        int centerY = currentObject.y + currentObject.height / 2;
        for (size_t k = 0; k < it->second.size(); k++) {
          RealObjectInformations pos = it->second.at(k);
          int posX = i % 2 == 0 ? pos.leftX : pos.rightX;
          int posY = i % 2 == 0 ? pos.leftY : pos.rightY;
          if (sqrt((centerX - posX) * (centerX - posX) + (centerY - posY) * (centerY - posY)) < 15) {
            this->squares.at(j).state = "1";
            break;
          }
        }
      }

      Utility::saveObjectsInformations(squarePath, this->squares);
    }
    Utility::drawSquares(p, this->squares, scale);
  }

  //ui->image->setCursor(Qt::UpArrowCursor);
  ui->image->setPixmap(rescaled);
}

void SecondaryObjectLocalisationWindow::on_pushButton_clicked()
{
  selectImage(selected + 1);
}

void SecondaryObjectLocalisationWindow::on_pushButton_2_clicked()
{
  selectImage(selected - 1);
}

bool SecondaryObjectLocalisationWindow::eventFilter(QObject *obj, QEvent *event)
{

  if (obj->objectName() == "image" && (event->type() == QEvent::MouseButtonRelease)) {
    QMouseEvent * mouseEvent = static_cast <QMouseEvent *> (event);

    if (mouseEvent->button() == 1) {
      if (this->localisationMode) {
        this->addSecondaryObject(mouseEvent->x(), mouseEvent->y());
      } else {
        this->switchSquare(mouseEvent->x(), mouseEvent->y());
      }
    } else {
      if (this->localisationMode) {
        this->removeSecondaryObject(mouseEvent->x(), mouseEvent->y());
      }
    }
  }

  if (event->type() == QEvent::KeyRelease) {
    QKeyEvent * keyEvent = static_cast <QKeyEvent *> (event);
    if (keyEvent->key() == Qt::Key_Left) {
      on_pushButton_2_clicked();
    }
    if (keyEvent->key() == Qt::Key_Right) {
      on_pushButton_clicked();
    }
  }
  return false;
}

void SecondaryObjectLocalisationWindow::addSecondaryObject(int x, int y) {
  SecondaryObjectPosition position;
  position.x = x / scale;
  position.y = y / scale;
  position.imageType = 1;

  secondaryPositions.push_back(position);
  saveSecondaryObjectsFile();
  selectImage(selected);
}

void SecondaryObjectLocalisationWindow::removeSecondaryObject(int x, int y) {
  x /= scale;
  y /= scale;

  float bestDistance = -1;
  int chosen = -1;

  for (size_t i = 0; i < secondaryPositions.size(); i++) {
    SecondaryObjectPosition pos = secondaryPositions.at(i);
    float distance = sqrt((pos.x - x) * (pos.x - x) + (pos.y - y) * (pos.y - y));
    if (chosen < 0 || distance < bestDistance) {
      bestDistance = distance;
      chosen = i;
    }
  }

  if (chosen != -1) {
    secondaryPositions.erase(secondaryPositions.begin() + chosen);
  }
  saveSecondaryObjectsFile();
  selectImage(selected);
}

void SecondaryObjectLocalisationWindow::saveSecondaryObjectsFile() {
  Utility::saveSecondaryObjectsFile(secondaryObjectsPath, secondaryPositions);
}

void SecondaryObjectLocalisationWindow::loadSecondaryObjectsFile() {
  secondaryPositions = Utility::loadSecondaryObjectsFile(secondaryObjectsPath);
}

void SecondaryObjectLocalisationWindow::switchSquare(int x, int y) {
  x /= scale;
  y /= scale;

  float bestDistance = -1;
  int chosen = -1;

  for (size_t i = 0; i < squares.size(); i++) {
    ObjectInformations square = squares.at(i);
    float distance = sqrt((square.x + square.width / 2 - x) * (square.x + square.width / 2 - x) + (square.y + square.height / 2 - y) * (square.y + square.height / 2 - y));
    if (chosen < 0 || distance < bestDistance) {
      bestDistance = distance;
      chosen = i;
    }
    if (distance > bestDistance - 0.001 && distance < bestDistance + 0.001 && rand() % 2 == 0) {
      chosen = i;
    }
  }

  if (chosen != -1) {
    squares.at(chosen).state = squares.at(chosen).state == "0" ? "1" : "0";
  }
  Utility::saveObjectsInformations(squarePath, this->squares);
  selectImage(selected);
}
