#ifndef SECONDARYOBJECTLOCALISATIONWINDOW_H
#define SECONDARYOBJECTLOCALISATIONWINDOW_H

#include <QMainWindow>
#include <QPixmap>
#include <QString>
#include <vector>
#include <map>
#include <QStringList>
#include "utility.h"
#include <QMouseEvent>
#include <QPainter>
#include <QByteArray>
#include <math.h>
#include "robonauteye.h"
#include "imageutility.h"

namespace Ui {
class SecondaryObjectLocalisationWindow;
}

class SecondaryObjectLocalisationWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit SecondaryObjectLocalisationWindow(QWidget *parent = 0, bool localisationMode = true);
  void selectImage(int i);
  ~SecondaryObjectLocalisationWindow();
  bool eventFilter(QObject *obj, QEvent *ev);
  void addSecondaryObject(int x, int y);
  void removeSecondaryObject(int x, int y);
  void switchSquare(int x, int y);
  void saveSecondaryObjectsFile();
  void loadSecondaryObjectsFile();

  
private slots:
  void on_pushButton_clicked();

  void on_pushButton_2_clicked();

private:
  Ui::SecondaryObjectLocalisationWindow *ui;
  map <QString, vector <RealObjectInformations> > infos;
  vector <QString> filesList;
  int selected;
  QString imagePath;
  QString secondaryObjectsPath;
  QString squarePath;
  vector <SecondaryObjectPosition> secondaryPositions;
  vector <ObjectInformations> squares;
  float scale;
  bool localisationMode;
};

#endif // SECONDARYOBJECTLOCALISATIONWINDOW_H
