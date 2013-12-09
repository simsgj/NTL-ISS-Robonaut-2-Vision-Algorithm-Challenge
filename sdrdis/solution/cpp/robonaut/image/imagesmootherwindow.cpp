#include "imagesmootherwindow.h"
#include "ui_imagesmootherwindow.h"

ImageSmootherWindow::ImageSmootherWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ImageSmootherWindow)
{
    ui->setupUi(this);

    this->sceneFrom = new QGraphicsScene(this);
    ui->viewFrom->setScene(this->sceneFrom);

    this->sceneTo = new QGraphicsScene(this);
    ui->viewTo->setScene(this->sceneTo);
}

ImageSmootherWindow::~ImageSmootherWindow()
{
    delete ui;
}

void ImageSmootherWindow::setImage(Image<int> * image) {
    this->image = image;
    this->refresh();
}

void ImageSmootherWindow::refresh() {
    float sigma = this->ui->sigmaSlider->value() / 100.0;
    this->ui->sigmaValue->setText(QString().setNum(sigma));

    Image<int> * dest = ImageSmoother<int>::apply(this->image, sigma);

    QImage imageFrom = ImageUtility::imageToQImage(this->image);
    QPixmap pixmapFrom = QPixmap::fromImage(imageFrom);

    QImage imageTo = ImageUtility::imageToQImage(dest);
    QPixmap pixmapTo = QPixmap::fromImage(imageTo);

    delete dest;

    this->sceneFrom->clear();
    this->sceneFrom->addPixmap(pixmapFrom);
    this->sceneFrom->setSceneRect(0, 0, pixmapFrom.width(), pixmapFrom.height());

    this->sceneTo->clear();
    this->sceneTo->addPixmap(pixmapTo);
    this->sceneTo->setSceneRect(0, 0, pixmapTo.width(), pixmapTo.height());
}

void ImageSmootherWindow::on_sigmaSlider_valueChanged(int value)
{
    this->refresh();
}
