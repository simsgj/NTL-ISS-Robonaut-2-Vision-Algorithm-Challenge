/********************************************************************************
** Form generated from reading UI file 'overallwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OVERALLWINDOW_H
#define UI_OVERALLWINDOW_H

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

class Ui_OverallWindow
{
public:
    QWidget *centralwidget;
    QLabel *leftImage;
    QLabel *rightImage;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *OverallWindow)
    {
        if (OverallWindow->objectName().isEmpty())
            OverallWindow->setObjectName(QStringLiteral("OverallWindow"));
        OverallWindow->resize(1600, 850);
        centralwidget = new QWidget(OverallWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        leftImage = new QLabel(centralwidget);
        leftImage->setObjectName(QStringLiteral("leftImage"));
        leftImage->setGeometry(QRect(0, 0, 800, 800));
        leftImage->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        rightImage = new QLabel(centralwidget);
        rightImage->setObjectName(QStringLiteral("rightImage"));
        rightImage->setGeometry(QRect(800, 0, 800, 800));
        rightImage->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        OverallWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(OverallWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1600, 22));
        OverallWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(OverallWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        OverallWindow->setStatusBar(statusbar);

        retranslateUi(OverallWindow);

        QMetaObject::connectSlotsByName(OverallWindow);
    } // setupUi

    void retranslateUi(QMainWindow *OverallWindow)
    {
        OverallWindow->setWindowTitle(QApplication::translate("OverallWindow", "MainWindow", 0));
        leftImage->setText(QApplication::translate("OverallWindow", "TextLabel", 0));
        rightImage->setText(QApplication::translate("OverallWindow", "TextLabel", 0));
    } // retranslateUi

};

namespace Ui {
    class OverallWindow: public Ui_OverallWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OVERALLWINDOW_H
