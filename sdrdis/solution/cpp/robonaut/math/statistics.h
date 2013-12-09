#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <cmath>

using namespace std;

template<class T>
class Statistics
{
public:
  Statistics();

  static float mean(vector <T> &list);
  static float mean(T * list, size_t nb);

  static float variance(vector <T> &list);
  static float variance(T * list, size_t nb);

  static void normalize(vector <T> &list, float newMean = 0, float newSigma = 1);
  static void normalize(T * list, size_t nb, float newMean = 0, float newSigma = 1);

  static void applyLowerBound(vector <T> &list, T lowerBound);
  static void applyLowerBound(T * list, size_t nb, T lowerBound);

  static void applyUpperBound(vector <T> &list, T lowerBound);
  static void applyUpperBound(T * list, size_t nb, T lowerBound);
};

template<class T>
Statistics<T>::Statistics()
{
}


template<class T>
float Statistics<T>::mean(vector <T> &list) {
  return mean(&list[0], list.size());
}

template<class T>
float Statistics<T>::mean(T * list, size_t nb) {
  float mean = 0;
  for (size_t i = 0; i < nb; i++) {
    mean += list[i] * 1.0 / nb;
  }
  return mean;
}

template<class T>
float Statistics<T>::variance(vector <T> &list) {
  return variance(&list[0], list.size());
}

template<class T>
float Statistics<T>::variance(T * list, size_t nb) {
  float meanList = mean(list, nb);
  float variance = 0;
  for (size_t i = 0; i < nb; i++) {
    float diff = list[i] - meanList;
    variance += (diff * diff) * 1.0 / (nb - 1);
  }
  return variance;
}

template<class T>
void Statistics<T>::normalize(vector <T> &list, float newMean, float newSigma) {
  normalize(&list[0], list.size(), newMean, newSigma);
}

template<class T>
void Statistics<T>::normalize(T * list, size_t nb, float newMean, float newSigma) {
  float meanList = mean(list, nb); // twice same process (variance need mean)
  float sigma = sqrt(variance(list, nb));

  for (size_t i = 0; i < nb; i++) {
    list[i] = ((list[i] - meanList) * newSigma / sigma) + newMean;
  }
}

template<class T>
void Statistics<T>::applyLowerBound(vector <T> &list, T lowerBound) {
  applyLowerBound(&list[0], list.size(), lowerBound);
}

template<class T>
void Statistics<T>::applyLowerBound(T * list, size_t nb, T lowerBound) {
  for (size_t i = 0; i < nb; i++) {
    if (list[i] < lowerBound) {
      list[i] = lowerBound;
    }
  }
}

template<class T>
void Statistics<T>::applyUpperBound(vector <T> &list, T upperBound) {
  applyUpperBound(&list[0], list.size(), upperBound);
}

template<class T>
void Statistics<T>::applyUpperBound(T * list, size_t nb, T upperBound) {
  for (size_t i = 0; i < nb; i++) {
    if (list[i] > upperBound) {
      list[i] = upperBound;
    }
  }
}

#endif // STATISTICS_H
