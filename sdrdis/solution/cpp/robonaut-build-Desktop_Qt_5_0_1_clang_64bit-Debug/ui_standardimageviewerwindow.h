/********************************************************************************
** Form generated from reading UI file 'standardimageviewerwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_STANDARDIMAGEVIEWERWINDOW_H
#define UI_STANDARDIMAGEVIEWERWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_StandardImageViewerWindow
{
public:
    QMenuBar *menubar;
    QWidget *centralwidget;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *StandardImageViewerWindow)
    {
        if (StandardImageViewerWindow->objectName().isEmpty())
            StandardImageViewerWindow->setObjectName(QStringLiteral("StandardImageViewerWindow"));
        StandardImageViewerWindow->resize(800, 600);
        menubar = new QMenuBar(StandardImageViewerWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        StandardImageViewerWindow->setMenuBar(menubar);
        centralwidget = new QWidget(StandardImageViewerWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        StandardImageViewerWindow->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(StandardImageViewerWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        StandardImageViewerWindow->setStatusBar(statusbar);

        retranslateUi(StandardImageViewerWindow);

        QMetaObject::connectSlotsByName(StandardImageViewerWindow);
    } // setupUi

    void retranslateUi(QMainWindow *StandardImageViewerWindow)
    {
        StandardImageViewerWindow->setWindowTitle(QApplication::translate("StandardImageViewerWindow", "MainWindow", 0));
    } // retranslateUi

};

namespace Ui {
    class StandardImageViewerWindow: public Ui_StandardImageViewerWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_STANDARDIMAGEVIEWERWINDOW_H
