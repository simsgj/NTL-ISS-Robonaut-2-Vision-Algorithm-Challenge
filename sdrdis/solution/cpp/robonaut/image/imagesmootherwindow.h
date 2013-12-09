#ifndef IMAGESMOOTHERWINDOW_H
#define IMAGESMOOTHERWINDOW_H

#include <QMainWindow>
#include "image.h"
#include <QGraphicsScene>
#include "imageutility.h"
#include "imagesmoother.h"
#include <QDebug>

namespace Ui {
class ImageSmootherWindow;
}

class ImageSmootherWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit ImageSmootherWindow(QWidget *parent = 0);
    ~ImageSmootherWindow();
    void setImage(Image<int> * image);
    void refresh();
    
private slots:
    void on_sigmaSlider_valueChanged(int value);


private:
    Ui::ImageSmootherWindow *ui;
    Image<int> * image;
    QGraphicsScene * sceneFrom;
    QGraphicsScene * sceneTo;
};

#endif // IMAGESMOOTHERWINDOW_H
