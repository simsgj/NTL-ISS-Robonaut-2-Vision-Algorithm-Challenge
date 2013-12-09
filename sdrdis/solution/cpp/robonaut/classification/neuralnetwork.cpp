#include "neuralnetwork.h"

#define sigmoid(value) (1.0 / (1.0 + exp(-value)))

NeuralNetwork::NeuralNetwork() {

}

float * NeuralNetwork::evaluate(float * inputLayer) {
    float * intermediary = inputLayer;
    float * oldIntermediary;
    for (unsigned int i = 1; i < nbLayers; i++) {
      oldIntermediary = intermediary;
      intermediary = evaluateLayer(i, oldIntermediary); // could be optimized
      if (i > 1) { // don't delete first layer (image)
        delete oldIntermediary;
      }
    }
    return intermediary;
}

float NeuralNetwork::evaluateSingle(float * inputLayer) {
  float * result = evaluate(inputLayer);

  float singleResult = result[1];
  delete result;
  return singleResult;
}

float * NeuralNetwork::evaluateLayer(int layer, float * previousLayer) {
    unsigned int nbPreviousUnits  = layers[layer - 1];
    unsigned int nbUnits          = layers[layer];
    float * values = new float[nbUnits];

    for (unsigned int j = 0; j < nbUnits; j++) {
        values[j] = 0;

        values[j] += this->weights[layer - 1][0 + j] * 1.0;
        for (unsigned int i = 0; i < nbPreviousUnits; i++) {
            values[j] += this->weights[layer - 1][(i + 1) * nbUnits + j] * previousLayer[i];
        }
        values[j] = sigmoid(values[j]);
    }
    return values;
}
