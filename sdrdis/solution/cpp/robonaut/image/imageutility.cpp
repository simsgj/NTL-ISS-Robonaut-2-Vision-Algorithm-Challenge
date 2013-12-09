#include "imageutility.h"

ImageUtility::ImageUtility()
{
}

Image<int> * ImageUtility::qImageToImage(QImage image) {
  int * convertedImage = new int[image.width() * image.height()];
  int W = image.width();
  int H = image.height();

  for (int y = 0; y < H; y++) {
    for (int x = 0; x < W; x++) {
      convertedImage[x + y * W] = ((qBlue(image.pixel(x, y)) + qRed(image.pixel(x, y)) + qGreen(image.pixel(x, y))) / 3);
    }
  }

  Image<int> * res = new Image<int>(H, W, convertedImage, 1, false);

  return res;
}

Image<float> * ImageUtility::qImageToImageFloat(QImage image) {
  float * convertedImage = new float[image.width() * image.height()];
  int W = image.width();
  int H = image.height();

  for (int y = 0; y < H; y++) {
    for (int x = 0; x < W; x++) {
      convertedImage[x + y * W] = ((qBlue(image.pixel(x, y)) + qRed(image.pixel(x, y)) + qGreen(image.pixel(x, y))) / 3);
    }
  }

  Image<float> * res = new Image<float>(H, W, convertedImage, 1, false);

  return res;
}

QImage ImageUtility::imageToQImage(Image<int> * image) {
  QImage img(image->getWidth(), image->getHeight(), QImage::Format_RGB32);
  for (int x = 0; x < image->getWidth(); x++) {
    for (int y = 0; y < image->getHeight(); y++) {
      int color = imgPixel(image, x, y);
      img.setPixel(x, y, qRgb(color, color, color));
    }
  }
  return img;
}

QImage ImageUtility::imageFloatToQImage(Image<float> * image) {
  QImage img(image->getWidth(), image->getHeight(), QImage::Format_RGB32);
  for (int x = 0; x < image->getWidth(); x++) {
    for (int y = 0; y < image->getHeight(); y++) {
      int color = imgPixel(image, x, y);
      img.setPixel(x, y, qRgb(color, color, color));
    }
  }
  return img;
}

QPixmap ImageUtility::imageToQPixmap(Image<int> * image) {
  QImage img = imageToQImage(image);
  return QPixmap::fromImage(img);
}

QPixmap ImageUtility::imageFloatToQPixmap(Image<float> * image) {
  QImage img = imageFloatToQImage(image);
  return QPixmap::fromImage(img);
}


template<class T>
const T& ImageUtility::kClamp( const T& x, const T& low, const T& high )
{
    if ( x < low )       return low;
    else if ( high < x ) return high;
    else                 return x;
}

inline
int ImageUtility::changeBrightness( int value, int brightness )
    {
    return kClamp( value + brightness * 255 / 100, 0, 255 );
    }

inline
int ImageUtility::changeContrast( int value, int contrast )
    {
    return kClamp((( value - 127 ) * contrast / 100 ) + 127, 0, 255 );
    }

inline
int ImageUtility::changeGamma( int value, int gamma )
    {
    return kClamp( int( pow( value / 255.0, 100.0 / gamma ) * 255 ), 0, 255 );
    }

inline
int ImageUtility::changeUsingTable( int value, const int table[] )
    {
    return table[ value ];
    }

template< int operation( int, int ) >
QImage ImageUtility::changeImage( const QImage& image, int value )
    {
    QImage im = image;
    im.detach();
    if( im.colorCount() == 0 ) /* truecolor */
        {
        if( im.format() != QImage::Format_RGB32 ) /* just in case */
            im = im.convertToFormat( QImage::Format_RGB32 );
        int table[ 256 ];
        for( int i = 0;
             i < 256;
             ++i )
            table[ i ] = operation( i, value );
        if( im.hasAlphaChannel() )
            {
            for( int y = 0;
                 y < im.height();
                 ++y )
                {
                QRgb* line = reinterpret_cast< QRgb* >( im.scanLine( y ));
                for( int x = 0;
                     x < im.width();
                     ++x )
                    line[ x ] = qRgba( changeUsingTable( qRed( line[ x ] ), table ),
                        changeUsingTable( qGreen( line[ x ] ), table ),
                        changeUsingTable( qBlue( line[ x ] ), table ),
                        changeUsingTable( qAlpha( line[ x ] ), table ));
                }
            }
        else
            {
            for( int y = 0;
                 y < im.height();
                 ++y )
                {
                QRgb* line = reinterpret_cast< QRgb* >( im.scanLine( y ));
                for( int x = 0;
                     x < im.width();
                     ++x )
                    line[ x ] = qRgb( changeUsingTable( qRed( line[ x ] ), table ),
                        changeUsingTable( qGreen( line[ x ] ), table ),
                        changeUsingTable( qBlue( line[ x ] ), table ));
                }
            }
        }
    else
        {
        QVector<QRgb> colors = im.colorTable();
        for( int i = 0;
             i < im.colorCount();
             ++i )
            colors[ i ] = qRgb( operation( qRed( colors[ i ] ), value ),
                operation( qGreen( colors[ i ] ), value ),
                operation( qBlue( colors[ i ] ), value ));
        }
    return im;
    }


// brightness is multiplied by 100 in order to avoid floating point numbers
QImage ImageUtility::changeBrightness( const QImage& image, int brightness )
    {
    if( brightness == 0 ) // no change
        return image;
    return changeImage< changeBrightness >( image, brightness );
    }


// contrast is multiplied by 100 in order to avoid floating point numbers
QImage ImageUtility::changeContrast( const QImage& image, int contrast )
    {
    if( contrast == 100 ) // no change
        return image;
    return changeImage< changeContrast >( image, contrast );
    }

// gamma is multiplied by 100 in order to avoid floating point numbers
QImage ImageUtility::changeGamma( const QImage& image, int gamma )
    {
    if( gamma == 100 ) // no change
        return image;
    return changeImage< changeGamma >( image, gamma );
    }
