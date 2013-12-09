#ifndef IMAGE_H
#define IMAGE_H


#include <vector>
#include <math.h>
//#include <QDebug>
//#include "point.h"
using namespace std;

#define imgPixel(im, x, y) (im->img[x + y * im->width])

template <class T>
class Image
{
public:
    Image(int H, int W);
    Image(int H, int W, T * img, double scale = 1, bool external = false);
    Image (const Image<T> &);
    ~Image();

    static Image<T> * resize(Image<T> * from, int size);
    static Image<T> * resize(Image<T> * from, int newH, int newW);
    static Image<T> * extract(Image<T> * image, int fromX, int fromY, int fromH, int fromW);
    static Image<T> * reshape(Image<T> * image, int fromX, int fromY, int fromH, int fromW, int toH, int toW);
    static Image<T> * reshapeInWhite(Image<T> * image, int fromX, int fromY, int fromH, int fromW, int toH, int toW, bool autoRotate = true);
    static vector <Image<T>  > * reshapeRotated(Image<T> * image, int fromX, int fromY, int fromH, int fromW, int toH, int toW);
    Image<T> * getWindow(int x, int y, int imgSize, int toSize);
    int getWidth();
    int getHeight();
    double getScale();
    T * getVector();
    void generateProfile();
    T pixel(int x, int y);
    void setPixel(int x, int y, T value);
    float diff(int xFrom, int yFrom, int xTo, int yTo);

    //vector <Point> * getPoints(int threshold = 170, double header = 0.1);


    double scale;
    int width;
    int height;
    T * img;

    bool external;
};







/* We moved out of cpp file since we are on template...
 * See: http://stackoverflow.com/a/10375372/888390
 */

template <class T>
Image<T>::Image(int H, int W)
{
  this->height = H;
  this->width = W;
  this->img = new T[this->height * this->width];
  this->scale = 1;
  this->external = false;
}


template <class T>
Image<T>::Image(int H, int W, T * img, double scale, bool external)
{
  this->height = H;
  this->width = W;
  this->img = img;
  this->scale = scale;
  this->external = external;
}

template <class T>
Image<T>::Image(const Image<T> &from)
{
  this->width = from.width;
  this->height = from.height;
  this->scale = from.scale;
  this->external = from.external;
  this->img = new T[this->width * this->height];
  size_t length = this->width * this->height;
  for (size_t i = 0; i < length; i++) {
    this->img[i] = from.img[i];
  }
}

template <class T>
Image<T>::~Image() {
  if (!this->external) {
    delete img;
  }
}

template <class T>
Image<T> * Image<T>::resize(Image<T> * from, int size) {
  double scale;
  if (from->height > from->width) {
    scale = size * 1.0 / from->getHeight();
  } else {
    scale = size * 1.0 / from->getWidth();
  }

  Image * img = Image<T>::resize(from, scale * from->getHeight(), scale * from->getWidth());
  img->scale = scale;
  return img;
}

template <class T>
Image<T> * Image<T>::resize(Image<T> * from, int newH, int newW) {
  T * newImage = new T[newH * newW];

  double rapX = from->width * 1.0 / newW;
  double rapY = from->height * 1.0 / newH;
  for (int y = 0; y < newH; y++) {
    for (int x = 0; x < newW; x++) {
      int newLeft = max((int)(x * rapX), 0);
      int newRight = min((int)(newLeft + rapX - 1), from->getWidth() - 1);
      int newTop = max((int)(y * rapY), 0);
      int newBottom = min((int)(newTop + rapY - 1), from->getHeight() - 1);

      T pixelValue = (imgPixel(from, newLeft, newTop) + imgPixel(from, newRight, newTop) + imgPixel(from, newLeft, newBottom) + imgPixel(from, newRight, newBottom)) / 4;
      //qDebug() << newRight + newBottom * W << image->size();
      newImage[x + y * newW] = pixelValue;
    }
  }
  return new Image<T>(newH, newW, newImage);
}

template <class T>
Image<T> * Image<T>::extract(Image<T> * image, int fromX, int fromY, int fromH, int fromW) {
  T * newImage = new T[fromH * fromW];

  for (int y = 0; y < fromH; y++) {
    for (int x = 0; x < fromW; x++) {
      newImage[x + y * fromW] = imgPixel(image, x + fromX, y + fromY);
    }
  }
  return new Image<T>(fromH, fromW, newImage);
}

template <class T>
Image<T> * Image<T>::reshape(Image<T> * image, int fromX, int fromY, int fromH, int fromW, int toH, int toW) {
  T * newImage = new T[toH * toW];

  double rapX = fromW * 1.0 / toW;
  double rapY = fromH * 1.0 / toH;
  for (int y = 0; y < toH; y++) {
    for (int x = 0; x < toW; x++) {
      int newLeft = max((int)(x * rapX + fromX), 0);
      int newRight = min((int)(newLeft + rapX - 1), fromX + fromW - 1);
      int newTop = max((int)(y * rapY + fromY), 0);
      int newBottom = min((int)(newTop + rapY - 1), fromY + fromH - 1);

      T pixelValue = (imgPixel(image, newLeft, newTop) + imgPixel(image, newRight, newTop) + imgPixel(image, newLeft, newBottom) + imgPixel(image, newRight, newBottom)) / 4;
      newImage[x + y * toW] = pixelValue; // can be highly optimized
    }
  }
  return new Image<T>(toH, toW, newImage);
}

template <class T>
Image<T> * Image<T>::reshapeInWhite(Image<T> * image, int fromX, int fromY, int fromH, int fromW, int toH, int toW, bool autoRotate) {
  size_t length = toH * toW;
  T * newImage = new T[length];

  for (size_t i = 0; i < length; i++) {
    newImage[i] = 120;
  }
  bool rotate = false;

  if (fromH > fromW && autoRotate) {
    rotate = true;
    int temp = toH;
    toH = toW;
    toW = temp;
  }

  double rapX = fromW * 1.0 / toW;
  double rapY = fromH * 1.0 / toH;


  double rap = max(rapX, rapY);
  int translationX = rapY > rapX ? ((fromW * 1.0) - toW * rap) / 2 : 0;
  int translationY = rapY > rapX ? 0 : ((fromH * 1.0) - toH * rap) / 2;

  for (int y = 0; y < toH; y++) {
    for (int x = 0; x < toW; x++) {
      int newLeft = (x * rap + fromX) + translationX;
      int newRight = (newLeft + rap - 1);
      int newTop = (y * rap + fromY) + translationY;
      int newBottom = (newTop + rap - 1);

      int leftTopPixel = (newLeft < fromX || newTop < fromY || newLeft > fromX + fromW - 1 || newTop > fromY + fromH - 1) ? 255 : imgPixel(image, newLeft, newTop);
      int rightTopPixel = (newRight < fromX || newTop < fromY || newRight > fromX + fromW - 1 || newTop > fromY + fromH - 1) ? 255 : imgPixel(image, newRight, newTop);
      int leftBottomPixel = (newLeft < fromX || newBottom < fromY || newLeft > fromX + fromW - 1 || newBottom > fromY + fromH - 1) ? 255 : imgPixel(image, newLeft, newBottom);
      int rightBottomPixel = (newRight < fromX || newBottom < fromY || newRight > fromX + fromW - 1 || newBottom > fromY + fromH - 1) ? 255 : imgPixel(image, newRight, newBottom);


      T pixelValue = (leftTopPixel + rightTopPixel + leftBottomPixel + rightBottomPixel) / 4;
      newImage[rotate ? (toH - y - 1 + x * toH) : (x + y * toW)] = pixelValue;
    }
  }
  return new Image<T>(rotate ? toW : toH, rotate ? toH : toW, newImage);
}

template <class T>
vector <Image <T> > * Image<T>::reshapeRotated(Image<T> * image, int fromX, int fromY, int fromH, int fromW, int toH, int toW) {
  vector <Image <T> > * images = new vector <Image <T> > ();
  size_t length = toH * toW;
  T * img0 = new T[length];
  T * img90 = new T[length];
  T * img180 = new T[length];
  T * img270 = new T[length];

  for (size_t i = 0; i < length; i++) {
    img0[i] = 120;
    img90[i] = 120;
    img180[i] = 120;
    img270[i] = 120;
  }


  double rapX = fromW * 1.0 / toW;
  double rapY = fromH * 1.0 / toH;


  double rap = max(rapX, rapY);
  int translationX = rapY > rapX ? ((fromW * 1.0) - toW * rap) / 2 : 0;
  int translationY = rapY > rapX ? 0 : ((fromH * 1.0) - toH * rap) / 2;

  for (int y = 0; y < toH; y++) {
    for (int x = 0; x < toW; x++) {
      int newLeft = (x * rap + fromX) + translationX;
      int newRight = (newLeft + rap - 1);
      int newTop = (y * rap + fromY) + translationY;
      int newBottom = (newTop + rap - 1);

      int leftTopPixel = (newLeft < fromX || newTop < fromY || newLeft > fromX + fromW - 1 || newTop > fromY + fromH - 1) ? 255 : imgPixel(image, newLeft, newTop);
      int rightTopPixel = (newRight < fromX || newTop < fromY || newRight > fromX + fromW - 1 || newTop > fromY + fromH - 1) ? 255 : imgPixel(image, newRight, newTop);
      int leftBottomPixel = (newLeft < fromX || newBottom < fromY || newLeft > fromX + fromW - 1 || newBottom > fromY + fromH - 1) ? 255 : imgPixel(image, newLeft, newBottom);
      int rightBottomPixel = (newRight < fromX || newBottom < fromY || newRight > fromX + fromW - 1 || newBottom > fromY + fromH - 1) ? 255 : imgPixel(image, newRight, newBottom);


      T pixelValue = (leftTopPixel + rightTopPixel + leftBottomPixel + rightBottomPixel) / 4;
      int rotX = toH - y - 1;
      int rotY = x;

      img0[x + y * toW] = pixelValue;
      img90[rotX + rotY * toH] = pixelValue;
      img180[toW - x - 1 + (toH - y - 1) * toW] = pixelValue;
      img270[toH - rotX - 1 + (toW - rotY - 1) * toH] = pixelValue;
    }
  }
  images->push_back(Image<T>(toH, toW, img0));
  images->push_back(Image<T>(toW, toH, img90));
  images->push_back(Image<T>(toH, toW, img180));
  images->push_back(Image<T>(toW, toH, img270));
  return images;
}

template <class T>
Image<T> * Image<T>::getWindow(int x, int y, int imgSize, int toSize) {
  /*
  int fromH;
  int fromW;
  vector < int > * img;

  img = this->img;
  imgSize = imgSize / sScale;
  fromH = this->H;
  fromW = this->W;
  */
  /*
  if (toSize > imgSize) {

  } else {
    img = this->sImg;
    fromH = this->sH;
    fromW = this->sW;
  }*/

  return reshape(this, x, y, imgSize, imgSize, toSize, toSize);
}

template <class T>
int Image<T>::getWidth() {
  return this->width;
}

template <class T>
int Image<T>::getHeight() {
  return this->height;
}

template <class T>
void Image<T>::generateProfile() {
  /*
  this->profileSum = 0;
  this->pH = this->height / 100;
  this->pW = this->width / 100;

  long threshold = PROFILE_THRESHOLD * PROFILE_SIZE * PROFILE_SIZE;


  this->profileXY.resize(this->sW + this->sH);
  this->pSum.resize(this->pW * this->pH, 0);
  this->pImg.resize(this->pW * this->pH, 0);


  for (long x = 0; x < this->sW; x++) {
    for (long y = 0; y < this->sH; y++) {

      long val = this->sImg->at(x + y * this->sW);

      this->profileXY[x] += val;

      this->profileXY[this->sW + y] += val;
      this->profileSum += val;


      long pIndex = (x / PROFILE_SIZE) + (y / PROFILE_SIZE) * this->pW;

      this->pSum[pIndex] += val;
      //qDebug() << pIndex << this->pSum[pIndex];
      if (this->pSum[pIndex] > threshold) {
        this->pImg[pIndex] = 255;
      }
    }
  }
  */
}

template <class T>
T Image<T>::pixel(int x, int y) {
  return this->img[min(x, this->width) + min(y, this->height) * this->width];
}

template <class T>
void Image<T>::setPixel(int x, int y, T value) {
  this->img[x + y * this->width] = value;
}

template <class T>
double Image<T>::getScale() {
  return this->scale;
}

template <class T>
T * Image<T>::getVector() {
  return this->img;
}

template <class T>
float Image<T>::diff(int xFrom, int yFrom, int xTo, int yTo) {
    float diff = imgPixel(this, xFrom, yFrom) - imgPixel(this, xTo, yTo);
    return sqrt(diff * diff + diff * diff + diff * diff);
}


#endif // IMAGE_H
