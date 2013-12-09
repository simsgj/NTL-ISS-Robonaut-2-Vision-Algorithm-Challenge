/********************************************************************************
** Form generated from reading UI file 'testprojectwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TESTPROJECTWINDOW_H
#define UI_TESTPROJECTWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_TestProjectWindow
{
public:
    QWidget *centralwidget;
    QLabel *image;
    QMenuBar *menubar;

    void setupUi(QMainWindow *TestProjectWindow)
    {
        if (TestProjectWindow->objectName().isEmpty())
            TestProjectWindow->setObjectName(QStringLiteral("TestProjectWindow"));
        TestProjectWindow->resize(1024, 800);
        centralwidget = new QWidget(TestProjectWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        image = new QLabel(centralwidget);
        image->setObjectName(QStringLiteral("image"));
        image->setGeometry(QRect(0, 0, 1024, 800));
        image->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        TestProjectWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(TestProjectWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1024, 22));
        TestProjectWindow->setMenuBar(menubar);

        retranslateUi(TestProjectWindow);

        QMetaObject::connectSlotsByName(TestProjectWindow);
    } // setupUi

    void retranslateUi(QMainWindow *TestProjectWindow)
    {
        TestProjectWindow->setWindowTitle(QApplication::translate("TestProjectWindow", "MainWindow", 0));
        image->setText(QApplication::translate("TestProjectWindow", "TextLabel", 0));
    } // retranslateUi

};

namespace Ui {
    class TestProjectWindow: public Ui_TestProjectWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TESTPROJECTWINDOW_H
