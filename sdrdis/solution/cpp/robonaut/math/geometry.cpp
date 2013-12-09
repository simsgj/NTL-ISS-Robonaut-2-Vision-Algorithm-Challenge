#include "geometry.h"

Geometry::Geometry()
{
}



OrthographicProjectionEquation Geometry::getOrthographicProjectionEquation(Point fromA, Point fromB, Point fromC, Point toA, Point toB, Point toC) {
  OrthographicProjectionEquation eq;
  eq.xEq = getOrthographicProjectionEquationItem(fromA, fromB, fromC, toA, toB, toC);
  swapPoint(fromA);
  swapPoint(fromB);
  swapPoint(fromC);
  swapPoint(toA);
  swapPoint(toB);
  swapPoint(toC);
  eq.yEq = getOrthographicProjectionEquationItem(fromA, fromB, fromC, toA, toB, toC);
  return eq;
}

OrthographicProjectionEquationItem Geometry::getOrthographicProjectionEquationItem(Point fromA, Point fromB, Point fromC, Point toA, Point toB, Point toC, bool tryToPreventFailing) {
  OrthographicProjectionEquationItem eq;
  float dividerB = ((2 * fromB.y - fromA.y - fromC.y) * (2 * fromA.x - fromB.x - fromC.x) - (2 * fromA.y - fromB.y - fromC.y) * (2 * fromB.x - fromA.x - fromC.x));
  float dividerA = (2 * fromA.x - fromB.x - fromC.x);

  if (dividerA == 0 || dividerB == 0) {
    if (tryToPreventFailing) {
      eq = getOrthographicProjectionEquationItem(fromA, fromC, fromB, toA, toC, toB, false);
      if (eq.correct) return eq;
      eq = getOrthographicProjectionEquationItem(fromB, fromA, fromC, toB, toA, toC, false);
      if (eq.correct) return eq;
      eq = getOrthographicProjectionEquationItem(fromB, fromC, fromA, toB, toC, toA, false);
      if (eq.correct) return eq;
      eq = getOrthographicProjectionEquationItem(fromC, fromA, fromB, toC, toA, toB, false);
      if (eq.correct) return eq;
      eq = getOrthographicProjectionEquationItem(fromC, fromB, fromA, toC, toB, toA, false);
      if (eq.correct) return eq;
    }
    eq.correct = false;
    return eq;
  }

  eq.correct = true;

  eq.b = ((2 * toB.x - toA.x - toC.x) * (2 * fromA.x - fromB.x - fromC.x) - (2 * toA.x - toB.x - toC.x) * (2 * fromB.x - fromA.x - fromC.x)) / dividerB;

  eq.a = (2 * toA.x - toB.x - toC.x - eq.b * (2 * fromA.y - fromB.y - fromC.y)) / dividerA;
  eq.c = (toA.x + toB.x + toC.x - eq.a * (fromA.x + fromB.x + fromC.x) - eq.b * (fromA.y + fromB.y + fromC.y)) / 3.0;
  return eq;
}



void Geometry::swapPoint(Point & from) {
  float tmp = from.x;
  from.x = from.y;
  from.y = tmp;
}

Point * Geometry::getOrthographicProjectionPoints(OrthographicProjectionEquation & eq, const Point * fromPoints, unsigned long nb) {
  Point * points = new Point[nb];
  for (unsigned long i = 0; i < nb; i++) {
    points[i].x = fromPoints[i].x * eq.xEq.a + fromPoints[i].y * eq.xEq.b + eq.xEq.c;
    points[i].y = fromPoints[i].y * eq.yEq.a + fromPoints[i].x * eq.yEq.b + eq.yEq.c;
  }
  return points;
}


// ([xyXY]\d)([(xyXY]) -> \1*\2
// (\))([xyXY]\d) -> \1*\2
// - -> -
// in fact I was searching for perspective projection
// http://stackoverflow.com/questions/8925569/perspective-projection-4-points
// http://mathematica.stackexchange.com/questions/9244/solve-system-of-equations-related-to-perspective-projection
// http://en.wikipedia.org/wiki/Projective_transformation
// http://mmlab.disi.unitn.it/wiki/index.php/2D_Homography:_Algorithm_and_Tool
// http://www.csse.uwa.edu.au/~pk/research/matlabfns/Projective/homography2d.m
// http://www.ics.forth.gr/~lourakis/homest/#Download
// http://stackoverflow.com/questions/3856072/svd-implementation-c
// http://math.nist.gov/tnt/download.html

// Test before if 4 points are possibles : distances between them compared to distance min / max, angle between them compared to distance min / max
// Point of secondary objects can be deduced from perspective projection + Point of points from models.csv file
// Other / complementary solution ; for 4 points, you make 4 sets of 3 points, and calculate orthographic calculation. a, b, and c should be similar (-+ epsilon), and contained within a min / max (as well as atan2(b, a) and b^2 + a^2 ?)
// -> One good point of this idea is that a part can be implemented before the 4th loop thus can prevent it (check if the 3 first point are correct)
// For collecting secondary object position matching the id, it should be possible to do a perspective projection on known points and get it. Then we can observe the calculated orthographic projections and know which limits we can set...
/*
 *
 *
 Perspective projection
 x' = (h00*x+h01*y+h02) / (h20*x+h21*y+1)
 y' = (h10*x+h11*y+h12) / (h20*x+h21*y+1)
 */
Array2D<double> * Geometry::getPerspectiveProjectionHomography(const Point * fromPoints, const Point * toPoints, unsigned long nb) {
  float fromScale = 0;
  float fromMeanX = 0;
  float fromMeanY = 0;
  float toScale = 0;
  float toMeanX = 0;
  float toMeanY = 0;
  const Point * fromPointsNormalized = fromPoints; //normalizePoints(fromPoints, nb, fromScale, fromMeanX, fromMeanY);
  const Point * toPointsNormalized = toPoints; //normalizePoints(toPoints, nb, toScale, toMeanX, toMeanY);

  Array2D<double> A(3 * nb, 9);

  for (unsigned long i = 0; i < nb; i++) {
    A[3 * i][0] = 0;
    A[3 * i][1] = 0;
    A[3 * i][2] = 0;
    A[3 * i][3] = -fromPointsNormalized[i].x;
    A[3 * i][4] = -fromPointsNormalized[i].y;
    A[3 * i][5] = -1;
    A[3 * i][6] = toPointsNormalized[i].y * fromPointsNormalized[i].x;
    A[3 * i][7] = toPointsNormalized[i].y * fromPointsNormalized[i].y;
    A[3 * i][8] = toPointsNormalized[i].y;

    A[3 * i + 1][0] = fromPointsNormalized[i].x;
    A[3 * i + 1][1] = fromPointsNormalized[i].y;
    A[3 * i + 1][2] = 1;
    A[3 * i + 1][3] = 0;
    A[3 * i + 1][4] = 0;
    A[3 * i + 1][5] = 0;
    A[3 * i + 1][6] = -toPointsNormalized[i].x * fromPointsNormalized[i].x;
    A[3 * i + 1][7] = -toPointsNormalized[i].x * fromPointsNormalized[i].y;
    A[3 * i + 1][8] = -toPointsNormalized[i].x;

    A[3 * i + 2][0] = -toPointsNormalized[i].y * fromPointsNormalized[i].x;
    A[3 * i + 2][1] = -toPointsNormalized[i].y * fromPointsNormalized[i].y;
    A[3 * i + 2][2] = -toPointsNormalized[i].y;
    A[3 * i + 2][3] = toPointsNormalized[i].x * fromPointsNormalized[i].x;
    A[3 * i + 2][4] = toPointsNormalized[i].x * fromPointsNormalized[i].y;
    A[3 * i + 2][5] = toPointsNormalized[i].x;
    A[3 * i + 2][6] = 0;
    A[3 * i + 2][7] = 0;
    A[3 * i + 2][8] = 0;
  }



  JAMA::SVD<double> svd(A);
  Array2D<double> VRaw;
  svd.getV(VRaw);

  Array2D<double> * H = new Array2D<double>(3, 3);
  (*H)[0][0] = VRaw[0][8];
  (*H)[0][1] = VRaw[1][8];
  (*H)[0][2] = VRaw[2][8];
  (*H)[1][0] = VRaw[3][8];
  (*H)[1][1] = VRaw[4][8];
  (*H)[1][2] = VRaw[5][8];
  (*H)[2][0] = VRaw[6][8];
  (*H)[2][1] = VRaw[7][8];
  (*H)[2][2] = VRaw[8][8];
/*
  Array2D<double> T1(3, 3);
  T1[0][0] = fromScale;
  T1[0][1] = 0;
  T1[0][2] = -fromMeanX * fromScale;
  T1[1][0] = 0;
  T1[1][1] = fromScale;
  T1[1][2] = -fromMeanY * fromScale;
  T1[2][0] = 0;
  T1[2][1] = 0;
  T1[2][2] = 1;

  Array2D<double> T2(3, 3);
  T2[0][0] = toScale;
  T2[0][1] = 0;
  T2[0][2] = -toMeanX * toScale;
  T2[1][0] = 0;
  T2[1][1] = toScale;
  T2[1][2] = -toMeanY * toScale;
  T2[2][0] = 0;
  T2[2][1] = 0;
  T2[2][2] = 1;

  Array2D<double> T2Inv(3, 3);
  T2Inv = matinvert(T2);
*/

  //(*H) = matmult(matmult(T2Inv, (*H)), T1);
  //qDebug() << "here" << (*H)[2][0];
  /*
  Point * res = new Point[nb];
  float meanMult = 0;
  for (unsigned long i = 0; i < nb; i++) {
    res[i].x = (*H)[0][0] * fromPoints[i].x + (*H)[0][1] * fromPoints[i].y + (*H)[0][2];
    res[i].y = (*H)[1][0] * fromPoints[i].x + (*H)[1][1] * fromPoints[i].y + (*H)[1][2];
    qDebug() << "rap" << toPoints[i].x / res[i].x;
    meanMult += toPoints[i].x / (res[i].x * nb * 2);
    meanMult += toPoints[i].y / (res[i].y * nb * 2);
  }


  for (unsigned long i = 0; i < 3; i++) {
    for (unsigned long j = 0; j < 3; j++) {
      (*H)[i][j] *= meanMult;
    }
  }
  */
/*
  delete fromPointsNormalized;
  delete toPointsNormalized;
*/
  return H;
}

Point * Geometry::normalizePoints(const Point * points, unsigned long nb, float & scale, float & meanX, float & meanY) {
  Point * newPoints = new Point[nb];
  meanX = 0;
  meanY = 0;
  for (unsigned long i = 0; i < nb; i++) {
    meanX += points[i].x / nb;
    meanY += points[i].y / nb;
  }
  float meanDistance = 0;
  for (unsigned long i = 0; i < nb; i++) {
    newPoints[i].x = points[i].x - meanX;
    newPoints[i].y = points[i].y - meanY;
    meanDistance += sqrt(newPoints[i].x * newPoints[i].x + newPoints[i].y * newPoints[i].y) / nb;
  }
  scale = sqrt(2) / meanDistance;
  for (unsigned long i = 0; i < nb; i++) {
    newPoints[i].x *= scale;
    newPoints[i].y *= scale;
  }
  return newPoints;
}

Point * Geometry::getPerspectiveProjectionPoints(Array2D<double> * H, const Point * fromPoints, unsigned long nb) {
  Point * toPoints = new Point[nb];
  for (unsigned long i = 0; i < nb; i++) {
    //qDebug() << "before" << (*H)[2][0];
    double w = (*H)[2][0] * fromPoints[i].x + (*H)[2][1] * fromPoints[i].y + (*H)[2][2];
    //qDebug() << w;
    toPoints[i].x = ((*H)[0][0] * fromPoints[i].x + (*H)[0][1] * fromPoints[i].y + (*H)[0][2]) / w;
    toPoints[i].y = ((*H)[1][0] * fromPoints[i].x + (*H)[1][1] * fromPoints[i].y + (*H)[1][2]) / w;
  }
  return toPoints;
}
