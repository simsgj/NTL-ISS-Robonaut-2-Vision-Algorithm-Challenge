#include "pca.h"

PCA::PCA()
{
}

float * PCA::reduce(float * inputLayer) {
    Statistics<float>::normalize(inputLayer, nbInputs);

    float * result = new float[nbEigenVectors];

    for (size_t i = 0; i < nbEigenVectors; i++) {
      result[i] = 0;
      for (size_t j = 0; j < nbInputs; j++) {
        result[i] += inputLayer[j] * this->eigenVectors[i * nbInputs + j];
      }
    }

    return result;
}
