#ifndef OVERALLWINDOW_H
#define OVERALLWINDOW_H

#include <QMainWindow>
#include <QEvent>
#include <QKeyEvent>
#include "utility.h"
#include "robonauteye.h"
#include <map>
#include <vector>
#include <math.h>

using namespace std;


namespace Ui {
class OverallWindow;
}

class OverallWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit OverallWindow(QWidget *parent = 0);
  ~OverallWindow();
  void loadInformations();
  void selectImage(int i);
  void drawFinalObjects(QPainter &leftPainter, QPainter &rightPainter, QString leftImagePath, QString rightImagePath, float scale);
  void drawPrimaryObjects(QPainter &p, QString imagePath, float scale);
  bool eventFilter(QObject *obj, QEvent *event);
  
private:
  Ui::OverallWindow *ui;
  map <QString, vector <RealObjectInformations> > infos;
  vector <QString> filesList;
  size_t selected;
};

#endif // OVERALLWINDOW_H
