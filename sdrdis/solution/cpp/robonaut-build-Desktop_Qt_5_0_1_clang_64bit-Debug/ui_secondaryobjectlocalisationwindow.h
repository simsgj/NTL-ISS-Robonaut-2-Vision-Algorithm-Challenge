/********************************************************************************
** Form generated from reading UI file 'secondaryobjectlocalisationwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SECONDARYOBJECTLOCALISATIONWINDOW_H
#define UI_SECONDARYOBJECTLOCALISATIONWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SecondaryObjectLocalisationWindow
{
public:
    QWidget *centralwidget;
    QLabel *image;
    QPushButton *pushButton;
    QPushButton *pushButton_2;

    void setupUi(QMainWindow *SecondaryObjectLocalisationWindow)
    {
        if (SecondaryObjectLocalisationWindow->objectName().isEmpty())
            SecondaryObjectLocalisationWindow->setObjectName(QStringLiteral("SecondaryObjectLocalisationWindow"));
        SecondaryObjectLocalisationWindow->resize(1024, 1007);
        centralwidget = new QWidget(SecondaryObjectLocalisationWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        image = new QLabel(centralwidget);
        image->setObjectName(QStringLiteral("image"));
        image->setGeometry(QRect(0, 0, 1024, 921));
        image->setCursor(QCursor(Qt::ArrowCursor));
        image->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        pushButton = new QPushButton(centralwidget);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setGeometry(QRect(910, 950, 114, 32));
        pushButton_2 = new QPushButton(centralwidget);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));
        pushButton_2->setGeometry(QRect(0, 950, 114, 32));
        SecondaryObjectLocalisationWindow->setCentralWidget(centralwidget);

        retranslateUi(SecondaryObjectLocalisationWindow);

        QMetaObject::connectSlotsByName(SecondaryObjectLocalisationWindow);
    } // setupUi

    void retranslateUi(QMainWindow *SecondaryObjectLocalisationWindow)
    {
        SecondaryObjectLocalisationWindow->setWindowTitle(QApplication::translate("SecondaryObjectLocalisationWindow", "MainWindow", 0));
        image->setText(QApplication::translate("SecondaryObjectLocalisationWindow", "Image", 0));
        pushButton->setText(QApplication::translate("SecondaryObjectLocalisationWindow", "Next", 0));
        pushButton_2->setText(QApplication::translate("SecondaryObjectLocalisationWindow", "Previous", 0));
    } // retranslateUi

};

namespace Ui {
    class SecondaryObjectLocalisationWindow: public Ui_SecondaryObjectLocalisationWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SECONDARYOBJECTLOCALISATIONWINDOW_H
