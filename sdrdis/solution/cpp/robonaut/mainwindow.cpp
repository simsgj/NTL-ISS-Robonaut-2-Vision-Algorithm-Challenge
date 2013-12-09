#include "mainwindow.h"
#include "ui_mainwindow.h"

vector <RealObjectInformations> MainWindow::secondaryPositions;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    secondaryObjectLocalisationWindow = new SecondaryObjectLocalisationWindow(this, true);
    secondaryObjectLocalisationSquaresWindow = new SecondaryObjectLocalisationWindow(this, false);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    QImage image(QString("data/test/grayscale.png"));
    imageWindow.setImage(ImageUtility::qImageToImage(image));
    imageWindow.show();
}

void MainWindow::on_pushButton_2_clicked()
{
    QImage image(QString("data/test/grayscale.png"));
    imageSmootherWindow.setImage(ImageUtility::qImageToImage(image));
    imageSmootherWindow.show();
}

void MainWindow::on_pushButton_3_clicked()
{
  std::vector<int> imageVector = Utility::loadImage("data/Lab2/LeftImage_ALJD.tif");
  // ISS/LeftImage_CIZA.tif
  // LOW contrast ISS/LeftImage_BHXG.tif
  // Sim/LeftImage_BGML.jpg
  // Lab/LeftImage_AHJH.tif
  // ISS/LeftImage_BHXG.tif
  // Lab2/LeftImage_ALJD.tif
  Image<float> * img = RobonautEye::convertVectorToImage(imageVector);
  Image<float> * imgResized = Image<float>::resize(img, 800);
  float ratio = imgResized->getWidth() * 1.0 / img->getWidth();
  delete img;

 map <QString, vector <RealObjectInformations> > infos = Utility::loadSecondaryObjects("data/model.csv", "data/secondary_positions/");
 secondaryPositions = infos.find("Lab2/ALJD.tif")->second;
 Utility::scaleInformations(secondaryPositions, ratio);


  segmentWindow.setImage(imgResized, MainWindow::isSetCorrect);
  segmentWindow.show();
}

void MainWindow::on_pushButton_4_clicked()
{
  secondaryObjectLocalisationWindow->selectImage(0);
  secondaryObjectLocalisationWindow->show();
}

bool MainWindow::isSetCorrect(Set<float> * set) {
  for (vector<RealObjectInformations>::iterator it = secondaryPositions.begin() ; it != secondaryPositions.end(); ++it) {
    float distance = sqrt(((*it).leftX - set->gravityCenterX) * ((*it).leftX - set->gravityCenterX) + ((*it).leftY - set->gravityCenterY) * ((*it).leftY - set->gravityCenterY));
    if (distance < 10) {
      return true;
    }
  }
  return false;
}

void MainWindow::on_pushButton_5_clicked()
{
  secondaryObjectLocalisationSquaresWindow->selectImage(0);
  secondaryObjectLocalisationSquaresWindow->show();
}

void MainWindow::on_pushButton_6_clicked()
{
  //testProjectWindow.setImage("data/ISS/", "CIZA.tif");
  testProjectWindow.show();
}

void MainWindow::on_pushButton_7_clicked()
{
  secondaryObjectClassifierWindow.selectImage(114);
  secondaryObjectClassifierWindow.show();
}

void MainWindow::on_pushButton_8_clicked()
{
    overallWindow.selectImage(57);
    overallWindow.show();
}
