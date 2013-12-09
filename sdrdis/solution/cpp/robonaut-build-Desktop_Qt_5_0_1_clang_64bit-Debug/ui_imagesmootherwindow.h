/********************************************************************************
** Form generated from reading UI file 'imagesmootherwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMAGESMOOTHERWINDOW_H
#define UI_IMAGESMOOTHERWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ImageSmootherWindow
{
public:
    QWidget *centralwidget;
    QGraphicsView *viewFrom;
    QGraphicsView *viewTo;
    QSlider *sigmaSlider;
    QLabel *label;
    QLabel *sigmaValue;
    QMenuBar *menubar;

    void setupUi(QMainWindow *ImageSmootherWindow)
    {
        if (ImageSmootherWindow->objectName().isEmpty())
            ImageSmootherWindow->setObjectName(QStringLiteral("ImageSmootherWindow"));
        ImageSmootherWindow->resize(630, 662);
        centralwidget = new QWidget(ImageSmootherWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        viewFrom = new QGraphicsView(centralwidget);
        viewFrom->setObjectName(QStringLiteral("viewFrom"));
        viewFrom->setGeometry(QRect(0, 0, 361, 321));
        viewTo = new QGraphicsView(centralwidget);
        viewTo->setObjectName(QStringLiteral("viewTo"));
        viewTo->setGeometry(QRect(0, 320, 361, 321));
        sigmaSlider = new QSlider(centralwidget);
        sigmaSlider->setObjectName(QStringLiteral("sigmaSlider"));
        sigmaSlider->setGeometry(QRect(420, 10, 160, 22));
        sigmaSlider->setMaximum(500);
        sigmaSlider->setValue(50);
        sigmaSlider->setOrientation(Qt::Horizontal);
        label = new QLabel(centralwidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(370, 10, 51, 21));
        sigmaValue = new QLabel(centralwidget);
        sigmaValue->setObjectName(QStringLiteral("sigmaValue"));
        sigmaValue->setGeometry(QRect(590, 10, 31, 21));
        ImageSmootherWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(ImageSmootherWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 630, 22));
        ImageSmootherWindow->setMenuBar(menubar);

        retranslateUi(ImageSmootherWindow);

        QMetaObject::connectSlotsByName(ImageSmootherWindow);
    } // setupUi

    void retranslateUi(QMainWindow *ImageSmootherWindow)
    {
        ImageSmootherWindow->setWindowTitle(QApplication::translate("ImageSmootherWindow", "MainWindow", 0));
        label->setText(QApplication::translate("ImageSmootherWindow", "Sigma:", 0));
        sigmaValue->setText(QApplication::translate("ImageSmootherWindow", "0", 0));
    } // retranslateUi

};

namespace Ui {
    class ImageSmootherWindow: public Ui_ImageSmootherWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMAGESMOOTHERWINDOW_H
