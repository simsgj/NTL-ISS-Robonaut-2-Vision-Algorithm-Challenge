#ifndef IMAGEWINDOW_H
#define IMAGEWINDOW_H

#include <QMainWindow>
#include "image.h"
#include "imageutility.h"
#include <QImage>
#include <QPixmap>
#include <QDebug>
#include <QGraphicsScene>


namespace Ui {
class ImageWindow;
}

class ImageWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit ImageWindow(QWidget *parent = 0);
    ~ImageWindow();
    void setImage(Image<int> * image);
    void refresh();
    
private slots:
    void on_sizeSlider_valueChanged(int value);

private:
    Ui::ImageWindow *ui;
    Image<int> * image;
    QGraphicsScene * scene;
};

#endif // IMAGEWINDOW_H
