#ifndef NEURALNETWORK_H
#define NEURALNETWORK_H

#include <vector>
#include <math.h>

using namespace std;

class NeuralNetwork
{
  public:
        NeuralNetwork();
        float * evaluate(float * inputLayer);
        float * evaluateLayer(int layer, float * previousLayer);
        float evaluateSingle(float * inputLayer);

        size_t * layers; // All layers
        size_t nbLayers;
        float ** weights;
        double precision;
        double recall;
        vector < vector <float> > probas;
};

#endif // NEURALNETWORK_H
