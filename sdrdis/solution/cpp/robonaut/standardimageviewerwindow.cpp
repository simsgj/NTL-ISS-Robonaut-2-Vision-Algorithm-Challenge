#include "standardimageviewerwindow.h"
#include "ui_standardimageviewerwindow.h"

StandardImageViewerWindow::StandardImageViewerWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::StandardImageViewerWindow)
{
  ui->setupUi(this);
  this->installEventFilter(this);
  loadInformations();
}

StandardImageViewerWindow::~StandardImageViewerWindow()
{
  delete ui;
}

void StandardImageViewerWindow::loadInformations() {
  this->infos = Utility::loadCSV("data/model.csv");
  for (map<QString, vector <RealObjectInformations> >::iterator it=infos.begin(); it!=infos.end(); ++it) {
    filesList.push_back(it->first);
  }
}

void StandardImageViewerWindow::selectImage(int i) {
  if (i < 0) {
    i = filesList.size() * 2 - 1;
  }
  if (i > (filesList.size() * 2 - 1)) { // after because unsigned
    i = 0;
  }
  selected = i;
  this->setWindowTitle(QString().setNum(i + 1) + " / " + QString().setNum(filesList.size() * 2));

  QStringList pathList = filesList.at(i / 2).split('/');
  pathList[1] = ((i % 2 == 0) ? "LeftImage_" : "RightImage_") + pathList[1];
  this->imagePath = "data/" + pathList.join('/');
}

bool StandardImageViewerWindow::eventFilter(QObject *obj, QEvent *event)
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

vector <Point> StandardImageViewerWindow::getPrimaryObjectsPositions() {
  vector <Point> positions;

  map <QString, vector <RealObjectInformations> >::iterator it = this->infos.find(filesList.at(selected / 2));
  positions = Utility::getPrimaryObjectsPositions(it->second, selected % 2 == 0);

  return positions;
}

vector <RealObjectInformations> StandardImageViewerWindow::getPrimaryObjectsRealInformations() {
  return this->infos.find(filesList.at(selected / 2))->second;
}
