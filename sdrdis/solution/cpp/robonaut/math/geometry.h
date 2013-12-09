#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "jama_svd.h"
#include "mathutility.h"
#include "tnt_extra.h"
using namespace std;

struct Point {
  double x;
  double y;

  Point() {
    this->x = 0;
    this->y = 0;
  }

  Point(float x, float y) {
    this->x = x;
    this->y = y;
  }
};

struct OrthographicProjectionEquationItem {
  float a;
  float b;
  float c;
  bool correct;
};

struct PerspectiveProjectionHomography {

};

struct OrthographicProjectionEquation {
  OrthographicProjectionEquationItem xEq;
  OrthographicProjectionEquationItem yEq;
};

class Geometry
{
public:
  Geometry();

  static OrthographicProjectionEquation getOrthographicProjectionEquation(Point fromA, Point fromB, Point fromC, Point toA, Point toB, Point toC);
  static void swapPoint(Point & from);
  static OrthographicProjectionEquationItem getOrthographicProjectionEquationItem(Point fromA, Point fromB, Point fromC, Point toA, Point toB, Point toC, bool tryToPreventFailing = true);
  static Point * getOrthographicProjectionPoints(OrthographicProjectionEquation & eq, const Point * fromPoints, unsigned long nb);

  static Array2D<double> * getPerspectiveProjectionHomography(const Point * fromPoints, const Point * toPoints, unsigned long nb);
  static Point * getPerspectiveProjectionPoints(Array2D<double> * H, const Point * fromPoints, unsigned long nb);
  static Point * normalizePoints(const Point * points, unsigned long nb, float & scale, float & meanX, float & meanY);
};

#endif // GEOMETRY_H
