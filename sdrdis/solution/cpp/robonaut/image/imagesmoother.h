#ifndef IMAGESMOOTHER_H
#define IMAGESMOOTHER_H

#define GAUSS_MASK_WIDTH 4.0

#include "image.h"
#include <cmath>

template <class T>
class ImageSmoother
{
public:
    ImageSmoother();

    static Image<T> * apply(Image<T> *src, float sigma);

    static vector<float> getGaussMask(float sigma);
    static void normalize(vector<float> &mask);
    static void convolveEven(Image<T> *src, Image<T> *res, vector<float> &mask);
};


template <class T>
ImageSmoother<T>::ImageSmoother()
{
}

template <class T>
vector<float> ImageSmoother<T>::getGaussMask(float sigma) {
    sigma = max(sigma, 0.01F);
    int maskLength = (int)ceil(sigma * GAUSS_MASK_WIDTH) + 1;
    vector<float> mask(maskLength);
    for (int i = 0; i < maskLength; i++) {
        mask[i] = exp(-0.5*i*i/(sigma*sigma));
    }
    return mask;
}

template <class T>
Image<T> * ImageSmoother<T>::apply(Image<T> *src, float sigma) {
  vector<float> mask = getGaussMask(sigma);
  normalize(mask);

  Image<T> *tmp = new Image<T>(src->getWidth(), src->getHeight());
  Image<T> *dst = new Image<T>(src->getHeight(), src->getWidth());
  convolveEven(src, tmp, mask);
  convolveEven(tmp, dst, mask);

  delete tmp;
  return dst;
}

template <class T>
void ImageSmoother<T>::normalize(vector<float> &mask) {
  int len = mask.size();
  float sum = 0;
  for (int i = 1; i < len; i++) {
    sum += fabs(mask[i]);
  }
  sum = 2*sum + fabs(mask[0]);
  for (int i = 0; i < len; i++) {
    mask[i] /= sum;
  }
}

template <class T>
void ImageSmoother<T>::convolveEven(Image<T> *src, Image<T> *res, vector<float> &mask) {
  int width = src->getWidth();
  int height = src->getHeight();
  int maskLength = mask.size();

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      float sum = mask[0] * imgPixel(src, x, y);
      for (int i = 1; i < maskLength; i++) {
        sum += mask[i] * (imgPixel(src, max(x-i,0), y) + imgPixel(src, min(x+i, width-1), y));
      }
      res->setPixel(y, x, sum);
    }
  }
}



#endif // IMAGESMOOTHER_H
