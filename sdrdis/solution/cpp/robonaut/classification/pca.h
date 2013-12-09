#ifndef PCA_H
#define PCA_H

#include <vector>
#include <math.h>
#include "statistics.h"

using namespace std;

class PCA
{
public:
    PCA();
    float * reduce(float * inputLayer);

protected:
    float * eigenVectors;
    size_t nbEigenVectors;
    size_t nbInputs;
};

#endif // PCA_H
