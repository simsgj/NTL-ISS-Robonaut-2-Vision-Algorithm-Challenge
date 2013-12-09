#ifndef TESTPROJECTWINDOW_H
#define TESTPROJECTWINDOW_H

#include <QMainWindow>
#include <QPixmap>
#include <vector>
#include <QEvent>
#include <QMouseEvent>
#include "standardimageviewerwindow.h"
#include "robonauteye.h"
#include "utility.h"
#include "geometry.h"
#include "tnt_array2d.h"

using namespace std;

namespace Ui {
class TestProjectWindow;
}

class TestProjectWindow : public StandardImageViewerWindow
{
  Q_OBJECT
  
public:
  explicit TestProjectWindow(QWidget *parent = 0);
  ~TestProjectWindow();
  void selectImage(int i);
  bool eventFilter(QObject *obj, QEvent *event);
  void addPoint(int x, int y);
  void removePoint(int x, int y);
  
private:
  Ui::TestProjectWindow *ui;
  float scale;
  vector <Point> positions;
  map <QString, vector <RealObjectInformations> > secondaryInfos;
};

#endif // TESTPROJECTWINDOW_H
