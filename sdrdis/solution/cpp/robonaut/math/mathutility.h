#ifndef MATHUTILITY_H
#define MATHUTILITY_H

#include <QDebug>
#include "tnt_array2d.h"

using namespace TNT;

class MathUtility
{
public:
  MathUtility();
  static void printArray2D(Array2D<double> & arr);
};

#endif // MATHUTILITY_H
