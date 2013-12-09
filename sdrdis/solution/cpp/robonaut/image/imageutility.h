#ifndef IMAGEUTILITY_H
#define IMAGEUTILITY_H

#include "image.h"
#include <vector>
#include <QImage>
#include <QPixmap>


using namespace std;
class ImageUtility
{
public:
    ImageUtility();

    static Image<int> * qImageToImage(QImage image);
    static Image<float> * qImageToImageFloat(QImage image);
    static QImage imageToQImage(Image<int> * image);
    static QPixmap imageToQPixmap(Image<int> * image);
    static QImage imageFloatToQImage(Image<float> * image);
    static QPixmap imageFloatToQPixmap(Image<float> * image);

    template<class T>
    static const T& kClamp( const T& x, const T& low, const T& high );

    static int changeBrightness( int value, int brightness );
    static int changeContrast( int value, int contrast );
    static int changeGamma( int value, int gamma );
    static int changeUsingTable( int value, const int table[] );

    template< int operation( int, int ) >
    static QImage changeImage( const QImage& image, int value );


    // brightness is multiplied by 100 in order to avoid floating point numbers
    static QImage changeBrightness( const QImage& image, int brightness );

    // contrast is multiplied by 100 in order to avoid floating point numbers
    static QImage changeContrast( const QImage& image, int contrast );

    // gamma is multiplied by 100 in order to avoid floating point numbers
    static QImage changeGamma( const QImage& image, int gamma );


};

#endif // IMAGEUTILITY_H
