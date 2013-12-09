#ifndef STANDARDIMAGEVIEWERWINDOW_H
#define STANDARDIMAGEVIEWERWINDOW_H

#include <QMainWindow>
#include <QEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <map>
#include <vector>
#include "utility.h"

using namespace std;

namespace Ui {
class StandardImageViewerWindow;
}

class StandardImageViewerWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit StandardImageViewerWindow(QWidget *parent = 0);
  ~StandardImageViewerWindow();
  void loadInformations();
  void selectImage(int i);
  vector <Point> getPrimaryObjectsPositions();
  vector <RealObjectInformations> getPrimaryObjectsRealInformations();
  bool eventFilter(QObject *obj, QEvent *ev);
  
protected:
  Ui::StandardImageViewerWindow *ui;
  map <QString, vector <RealObjectInformations> > infos;
  vector <QString> filesList;
  int selected;
  QString imagePath;
};

#endif // STANDARDIMAGEVIEWERWINDOW_H
