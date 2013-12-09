#include "mainwindow.h"
#include <QApplication>

#include "jama_svd.h"
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "tnt_math_utils.h"

#include <QElapsedTimer>
#include <QDebug>
#include "robonauteye.h"
#include "vector"
#include "utility.h"
#include "statistics.h"
#include "secondaryobjectlocalisationwindow.h"
#include <QImage>
#include "math.h"
#include "image.h"
#include "secondaryobjectpca.h"
#include "secondaryobjectneuralnetwork.h"
#include "geometry.h"
#include "tests.h"


using namespace std;

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  //Utility::getMostFrequentState("data/model.csv");
  /*Utility::exportDecals("data/decals");
  exit(0);*/
  /*testSwitchDecal();
  exit(0);*/
  //Utility::exportDecals("data/decals");
  /*Utility::exportPrimaryObjectsToImages("data/primary_object_classifier");
  exit(0);*/
  /*testSecondaryObjectPCAandNN();
  exit(0);*/
    /*Utility::getScoreStatistics();
    exit(0);*/
    /*Utility::exportObjectsInformationsToImages("data/squares", "data/squares_images");
    exit(0);*/
    // Utility::upgradeSecondaryPositions();
    /*
    vector<int> leftImage = Utility::loadImage("data/ISS/LeftImage_CIZA.tif");
    vector<int> rightImage = Utility::loadImage("data/ISS/RightImage_CIZA.tif");
    RobonautEye eye;
    vector<string> recognized = eye.recognizeObjects(leftImage, rightImage);
    for (size_t i = 0; i < recognized.size(); i++) {
      qDebug() << recognized[i].c_str();
    }
    exit(0);
    */
    MainWindow w;
    w.show();
    
    return a.exec();
}



