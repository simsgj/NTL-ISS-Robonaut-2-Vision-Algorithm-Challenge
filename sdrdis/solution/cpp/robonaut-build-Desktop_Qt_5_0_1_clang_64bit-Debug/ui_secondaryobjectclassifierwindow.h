/********************************************************************************
** Form generated from reading UI file 'secondaryobjectclassifierwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SECONDARYOBJECTCLASSIFIERWINDOW_H
#define UI_SECONDARYOBJECTCLASSIFIERWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SecondaryObjectClassifierWindow
{
public:
    QWidget *centralwidget;
    QLabel *image;
    QMenuBar *menubar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *SecondaryObjectClassifierWindow)
    {
        if (SecondaryObjectClassifierWindow->objectName().isEmpty())
            SecondaryObjectClassifierWindow->setObjectName(QStringLiteral("SecondaryObjectClassifierWindow"));
        SecondaryObjectClassifierWindow->resize(1024, 900);
        centralwidget = new QWidget(SecondaryObjectClassifierWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        image = new QLabel(centralwidget);
        image->setObjectName(QStringLiteral("image"));
        image->setGeometry(QRect(0, -20, 1024, 900));
        image->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        SecondaryObjectClassifierWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(SecondaryObjectClassifierWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1024, 22));
        SecondaryObjectClassifierWindow->setMenuBar(menubar);
        statusBar = new QStatusBar(SecondaryObjectClassifierWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        SecondaryObjectClassifierWindow->setStatusBar(statusBar);

        retranslateUi(SecondaryObjectClassifierWindow);

        QMetaObject::connectSlotsByName(SecondaryObjectClassifierWindow);
    } // setupUi

    void retranslateUi(QMainWindow *SecondaryObjectClassifierWindow)
    {
        SecondaryObjectClassifierWindow->setWindowTitle(QApplication::translate("SecondaryObjectClassifierWindow", "MainWindow", 0));
        image->setText(QApplication::translate("SecondaryObjectClassifierWindow", "image", 0));
    } // retranslateUi

};

namespace Ui {
    class SecondaryObjectClassifierWindow: public Ui_SecondaryObjectClassifierWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SECONDARYOBJECTCLASSIFIERWINDOW_H
