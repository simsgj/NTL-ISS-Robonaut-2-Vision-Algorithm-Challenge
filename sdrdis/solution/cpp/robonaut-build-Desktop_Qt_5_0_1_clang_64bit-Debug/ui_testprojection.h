/********************************************************************************
** Form generated from reading UI file 'testprojection.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TESTPROJECTION_H
#define UI_TESTPROJECTION_H

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

class Ui_testProjection
{
public:
    QWidget *centralwidget;
    QLabel *label;
    QLabel *label_2;
    QMenuBar *menubar;

    void setupUi(QMainWindow *testProjection)
    {
        if (testProjection->objectName().isEmpty())
            testProjection->setObjectName(QStringLiteral("testProjection"));
        testProjection->resize(1024, 800);
        centralwidget = new QWidget(testProjection);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        label = new QLabel(centralwidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(0, 0, 512, 512));
        label->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        label_2 = new QLabel(centralwidget);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(512, 0, 512, 512));
        label_2->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        testProjection->setCentralWidget(centralwidget);
        menubar = new QMenuBar(testProjection);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1024, 22));
        testProjection->setMenuBar(menubar);

        retranslateUi(testProjection);

        QMetaObject::connectSlotsByName(testProjection);
    } // setupUi

    void retranslateUi(QMainWindow *testProjection)
    {
        testProjection->setWindowTitle(QApplication::translate("testProjection", "MainWindow", 0));
        label->setText(QApplication::translate("testProjection", "TextLabel", 0));
        label_2->setText(QApplication::translate("testProjection", "TextLabel", 0));
    } // retranslateUi

};

namespace Ui {
    class testProjection: public Ui_testProjection {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TESTPROJECTION_H
