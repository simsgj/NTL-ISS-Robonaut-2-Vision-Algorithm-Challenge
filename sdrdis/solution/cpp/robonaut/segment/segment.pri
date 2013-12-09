#-------------------------------------------------
#
# Project created by QtCreator 2013-03-23T16:25:07
#
#-------------------------------------------------

INCLUDEPATH += $$PWD

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

HEADERS += \
    $$PWD/felzenhuttensegmentation.h \
    segment/segmentwindow.h \
    segment/imagedisjointset.h \
    segment/standardmerger.h \
    segment/segmentutility.h \
    segment/standardselector.h \
    segment/mergepair.h

SOURCES += \
    segment/segmentwindow.cpp \
    segment/segmentutility.cpp \
    segment/mergepair.cpp

FORMS += \
    segment/segmentwindow.ui
