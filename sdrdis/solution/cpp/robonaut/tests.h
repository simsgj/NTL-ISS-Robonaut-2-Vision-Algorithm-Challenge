#ifndef TESTS_H
#define TESTS_H

#include <set>
#include <QDebug>
#include "robonauteye.h"

using namespace std;

void testSwitchDecal() {
  Point decal = RobonautEye::getSwitchDecal(1, 1, 0, 3.1415 / 2 + 0.3);
  qDebug() << decal.x;
  qDebug() << decal.y;
}

void testGetQuadruples() {
  vector<ObjectInformations> points;
  points.push_back(ObjectInformations(0, 0));
  points.push_back(ObjectInformations(0, 5));
  points.push_back(ObjectInformations(4, 7));
  points.push_back(ObjectInformations(2, 1));
  points.push_back(ObjectInformations(6, 2));

  multiset<QuadrupleIndexes> * quadruples = RobonautEye::getQuadruples(points);
  multiset<QuadrupleIndexes>::iterator it;
  for (it=quadruples->begin(); it!=quadruples->end(); ++it) {
    qDebug() << (*it).i << (*it).j << (*it).k << (*it).l << (*it).score;
  }

  delete quadruples;
}

void testSecondaryObjectPCAandNN() {
  QImage imgRaw("data/squares_images/negatives/0.png");
  Image<float> * img = ImageUtility::qImageToImageFloat(imgRaw);

  SecondaryObjectPCA pca;
  float * reduced = pca.reduce(img->img);
  for (size_t i = 0; i < 20; i++) {
    qDebug() << i << img->img[i];
  }

  qDebug() << img->getWidth() << img->getHeight();
  for (size_t i = 0; i < 80; i++) {
    qDebug() << reduced[i];
  }

  SecondaryObjectNeuralNetwork nn;
  float * evaluate = nn.evaluate(reduced);
  qDebug() << "RES: " << evaluate[0] << evaluate[1] << evaluate[2] << evaluate[3];
  delete evaluate;
}

/*
    Point * pointsFrom = new Point[4];
    Point * pointsTo = new Point[4];

    pointsFrom[0].x = 100;
    pointsFrom[0].y = 2;
    pointsFrom[1].x = 3;
    pointsFrom[1].y = 4;
    pointsFrom[2].x = 5;
    pointsFrom[2].y = 6;
    pointsFrom[3].x = 20;
    pointsFrom[3].y = 8;

    pointsTo[0].x = 20;
    pointsTo[0].y = 30;
    pointsTo[1].x = 40;
    pointsTo[1].y = 50;
    pointsTo[2].x = 60;
    pointsTo[2].y = 70;
    pointsTo[3].x = 80;
    pointsTo[3].y = 90;

    QElapsedTimer timer;
    timer.start();
    Array2D<double> * H = Geometry::getPerspectiveProjectionHomography(pointsFrom, pointsTo, 4);

    qDebug() << "from";
    Point * pointsDeduced = Geometry::getPerspectiveProjectionPoints(H, pointsFrom, 4);
    for (size_t i = 0; i < 4; i++) {
      qDebug() << pointsTo[i].x << pointsTo[i].y << pointsDeduced[i].x << pointsDeduced[i].y;
    }

    return 0;
*/
/*
    Array2D<float> arrTest(2, 2);
    arrTest[0][0] = 1;
    arrTest[0][1] = 2;
    arrTest[1][0] = 3;
    arrTest[1][1] = 4;

    JAMA::SVD<float> svd(arrTest);
    Array2D<float> vTest;
    svd.getV(vTest);

    qDebug() << vTest[0][0] << vTest[0][1];
    qDebug() << vTest[1][0] << vTest[1][1];

    QElapsedTimer timer;
    timer.start();
    for (int x = 0; x < 1000; x++) {
        Array2D<float> arr(30, 30);
        for (int i = 0; i < 30; i++) {
            for (int j = 0; j < 30; j++) {
                arr[i][j] = i + j;
            }
        }

        JAMA::SVD<float> svd(arr);
        Array2D<float> v;
        svd.getV(v);
    }
    qDebug() << timer.elapsed() << "ms";
    exit(0);
*/



/*
    QImage imgRaw("data/squares_images/negatives/0.png");
    Image<float> * img = ImageUtility::qImageToImageFloat(imgRaw);
    SecondaryObjectPCA pca;
    float * reduced = pca.reduce(img->img);

    for (size_t i = 0; i < 20; i++) {
      qDebug() << i << img->img[i];
    }

    qDebug() << img->getWidth() << img->getHeight();
    for (size_t i = 0; i < 40; i++) {
      qDebug() << reduced[i];
    }

    SecondaryObjectNeuralNetwork nn;
    qDebug() << "RES: " << nn.evaluateSingle(reduced);

    exit(0);
*/

#endif // TESTS_H
