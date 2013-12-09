#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
#include "imagewindow.h"
#include "image.h"
#include "imageutility.h"
#include "imagesmootherwindow.h"
#include "segmentwindow.h"
#include "utility.h"
#include "robonauteye.h"
#include "secondaryobjectlocalisationwindow.h"
#include "felzenhuttensegmentation.h"
#include "testprojectwindow.h"
#include "secondaryobjectclassifierwindow.h"
#include "overallwindow.h"

using namespace std;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

  static bool isSetCorrect(Set<float> * set);
    
private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();


    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_8_clicked();

private:
    Ui::MainWindow *ui;
    ImageWindow imageWindow;
    ImageSmootherWindow imageSmootherWindow;
    SegmentWindow segmentWindow;
    SecondaryObjectLocalisationWindow * secondaryObjectLocalisationWindow;
    SecondaryObjectLocalisationWindow * secondaryObjectLocalisationSquaresWindow;
    TestProjectWindow testProjectWindow;
    SecondaryObjectClassifierWindow secondaryObjectClassifierWindow;
    OverallWindow overallWindow;
    static vector <RealObjectInformations> secondaryPositions;
};

#endif // MAINWINDOW_H
