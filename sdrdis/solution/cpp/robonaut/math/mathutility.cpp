#include "mathutility.h"

MathUtility::MathUtility()
{
}

void MathUtility::printArray2D(Array2D<double> & arr) {
  int sizeI = arr.dim1();
  int sizeJ = arr.dim2();
  qDebug() << "SIZE:" << sizeI << sizeJ;
  for (int i = 0; i < sizeI; i++) {
    QString str = "";
    for (int j = 0; j < sizeJ; j++) {
      str += QString().setNum((double)arr[i][j]) + ", ";
    }
    qDebug() << QString("[" + str + "]").toStdString().c_str();
  }
}
