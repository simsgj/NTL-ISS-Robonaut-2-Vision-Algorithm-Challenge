/********************************************************************************
** Form generated from reading UI file 'imagewindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMAGEWINDOW_H
#define UI_IMAGEWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QSlider>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ImageWindow
{
public:
    QWidget *centralWidget;
    QGraphicsView *graphicsView;
    QSlider *sizeSlider;
    QLabel *sizeLabel;
    QLabel *sizeValue;
    QToolBar *mainToolBar;

    void setupUi(QMainWindow *ImageWindow)
    {
        if (ImageWindow->objectName().isEmpty())
            ImageWindow->setObjectName(QStringLiteral("ImageWindow"));
        ImageWindow->resize(914, 523);
        centralWidget = new QWidget(ImageWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        graphicsView = new QGraphicsView(centralWidget);
        graphicsView->setObjectName(QStringLiteral("graphicsView"));
        graphicsView->setGeometry(QRect(0, 0, 551, 511));
        sizeSlider = new QSlider(centralWidget);
        sizeSlider->setObjectName(QStringLiteral("sizeSlider"));
        sizeSlider->setGeometry(QRect(620, 10, 241, 22));
        sizeSlider->setMinimum(0);
        sizeSlider->setMaximum(200);
        sizeSlider->setValue(100);
        sizeSlider->setOrientation(Qt::Horizontal);
        sizeLabel = new QLabel(centralWidget);
        sizeLabel->setObjectName(QStringLiteral("sizeLabel"));
        sizeLabel->setGeometry(QRect(560, 10, 62, 16));
        sizeValue = new QLabel(centralWidget);
        sizeValue->setObjectName(QStringLiteral("sizeValue"));
        sizeValue->setGeometry(QRect(870, 10, 41, 20));
        ImageWindow->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(ImageWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        ImageWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);

        retranslateUi(ImageWindow);

        QMetaObject::connectSlotsByName(ImageWindow);
    } // setupUi

    void retranslateUi(QMainWindow *ImageWindow)
    {
        ImageWindow->setWindowTitle(QApplication::translate("ImageWindow", "ImageWindow", 0));
        sizeLabel->setText(QApplication::translate("ImageWindow", "Size:", 0));
        sizeValue->setText(QApplication::translate("ImageWindow", "1", 0));
    } // retranslateUi

};

namespace Ui {
    class ImageWindow: public Ui_ImageWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMAGEWINDOW_H
