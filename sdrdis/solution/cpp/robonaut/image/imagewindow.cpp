#include "imagewindow.h"
#include "ui_imagewindow.h"

ImageWindow::ImageWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ImageWindow)
{
    ui->setupUi(this);
    this->scene = new QGraphicsScene(this);
    ui->graphicsView->setScene(this->scene);
}

ImageWindow::~ImageWindow()
{
    delete ui;
    delete this->scene;
}

void ImageWindow::setImage(Image<int> * image) {
    this->image = image;
    this->refresh();
}

void ImageWindow::refresh() {
    float ratio = this->ui->sizeSlider->value() * 2.0 / this->ui->sizeSlider->maximum();
    Image<int> * res = Image<int>::resize(this->image, this->image->height * ratio, this->image->width * ratio);
    QImage image = ImageUtility::imageToQImage(res);
    QPixmap pixmap = QPixmap::fromImage(image);


    this->scene->clear();
    this->scene->addPixmap(pixmap);
    this->scene->setSceneRect(0, 0, pixmap.width(), pixmap.height());

    this->ui->sizeValue->setText(QString().setNum(ratio));
    delete res;
}

void ImageWindow::on_sizeSlider_valueChanged(int value)
{
    refresh();
}
