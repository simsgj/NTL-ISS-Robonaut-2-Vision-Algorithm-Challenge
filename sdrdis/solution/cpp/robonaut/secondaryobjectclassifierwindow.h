#ifndef SECONDARYOBJECTCLASSIFIERWINDOW_H
#define SECONDARYOBJECTCLASSIFIERWINDOW_H

#include <QMainWindow>
#include <map>
#include <vector>
#include "robonauteye.h"
#include "utility.h"
#include <QEvent>
#include <QKeyEvent>
#include <QElapsedTimer>

using namespace std;

namespace Ui {
class SecondaryObjectClassifierWindow;
}

class SecondaryObjectClassifierWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit SecondaryObjectClassifierWindow(QWidget *parent = 0);
  ~SecondaryObjectClassifierWindow();
  void selectImage(int i);
  bool eventFilter(QObject *obj, QEvent *ev);
  
protected:
  Ui::SecondaryObjectClassifierWindow *ui;

  map <QString, vector <RealObjectInformations> > infos;
  vector <QString> filesList;
  int selected;
  float scale;
};

#endif // SECONDARYOBJECTCLASSIFIERWINDOW_H
