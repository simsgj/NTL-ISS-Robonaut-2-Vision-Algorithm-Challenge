#-------------------------------------------------
#
# Project created by QtCreator 2013-03-23T16:25:07
#
#-------------------------------------------------

INCLUDEPATH += $$PWD

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets


SOURCES += $$PWD/imagewindow.cpp \
    $$PWD/imageutility.cpp \
    $$PWD/imagesmootherwindow.cpp \
    image/imagesmoother.cpp

HEADERS  += $$PWD/imagewindow.h \
    $$PWD/image.h \
    $$PWD/imageutility.h \
    $$PWD/imagesmoother.h \
    $$PWD/imagesmootherwindow.h

FORMS    += $$PWD/imagewindow.ui \
    image/imagesmootherwindow.ui
