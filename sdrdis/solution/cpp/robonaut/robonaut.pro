#-------------------------------------------------
#
# Project created by QtCreator 2013-03-23T16:29:27
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = robonaut
TEMPLATE = app

include(math/math.pri)

SOURCES += main.cpp\
        mainwindow.cpp \
    robonauteye.cpp \
    utility.cpp \
    secondaryobjectlocalisationwindow.cpp \
    testprojectwindow.cpp \
    secondaryobjectpca.cpp \
    secondaryobjectneuralnetwork.cpp \
    secondaryobjectclassifierwindow.cpp \
    standardimageviewerwindow.cpp \
    primaryobjectpca.cpp \
    primaryobjectneuralnetwork.cpp \
    quadrupleindexes.cpp \
    overallwindow.cpp

HEADERS  += mainwindow.h \
    robonauteye.h \
    utility.h \
    secondaryobjectlocalisationwindow.h \
    testprojectwindow.h \
    secondaryobjectpca.h \
    secondaryobjectneuralnetwork.h \
    secondaryobjectclassifierwindow.h \
    standardimageviewerwindow.h \
    primaryobjectpca.h \
    primaryobjectneuralnetwork.h \
    quadrupleindexes.h \
    tests.h \
    overallwindow.h

FORMS    += mainwindow.ui \
    secondaryobjectlocalisationwindow.ui \
    testprojectwindow.ui \
    secondaryobjectclassifierwindow.ui \
    standardimageviewerwindow.ui \
    overallwindow.ui

include(image/image.pri)
include(graph/graph.pri)
include(segment/segment.pri)
include(classification/classification.pri)

OTHER_FILES += \
    garbage.txt
