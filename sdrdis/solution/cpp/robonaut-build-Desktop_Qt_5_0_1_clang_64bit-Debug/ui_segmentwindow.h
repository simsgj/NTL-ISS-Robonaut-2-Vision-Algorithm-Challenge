/********************************************************************************
** Form generated from reading UI file 'segmentwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SEGMENTWINDOW_H
#define UI_SEGMENTWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SegmentWindow
{
public:
    QWidget *centralwidget;
    QGraphicsView *viewFrom;
    QGraphicsView *viewSegmented;
    QGraphicsView *viewSelected1;
    QGraphicsView *viewMerged;
    QScrollArea *scrollArea;
    QWidget *scrollAreaWidgetContents;
    QGroupBox *groupBox_3;
    QWidget *verticalLayoutWidget_3;
    QVBoxLayout *selectLayout1;
    QGroupBox *groupBox;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *preprocessingLayout;
    QGroupBox *groupBox_2;
    QWidget *verticalLayoutWidget_2;
    QVBoxLayout *splitLayout;
    QGroupBox *groupBox_4;
    QWidget *verticalLayoutWidget_4;
    QVBoxLayout *mergeLayout;
    QGroupBox *groupBox_5;
    QWidget *verticalLayoutWidget_5;
    QVBoxLayout *selectLayout2;
    QGraphicsView *viewSelected2;
    QGraphicsView *viewRectangles;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *SegmentWindow)
    {
        if (SegmentWindow->objectName().isEmpty())
            SegmentWindow->setObjectName(QStringLiteral("SegmentWindow"));
        SegmentWindow->resize(1037, 982);
        centralwidget = new QWidget(SegmentWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        viewFrom = new QGraphicsView(centralwidget);
        viewFrom->setObjectName(QStringLiteral("viewFrom"));
        viewFrom->setGeometry(QRect(0, 0, 361, 321));
        viewSegmented = new QGraphicsView(centralwidget);
        viewSegmented->setObjectName(QStringLiteral("viewSegmented"));
        viewSegmented->setGeometry(QRect(360, 0, 361, 321));
        viewSelected1 = new QGraphicsView(centralwidget);
        viewSelected1->setObjectName(QStringLiteral("viewSelected1"));
        viewSelected1->setGeometry(QRect(0, 320, 361, 321));
        viewMerged = new QGraphicsView(centralwidget);
        viewMerged->setObjectName(QStringLiteral("viewMerged"));
        viewMerged->setGeometry(QRect(360, 320, 361, 321));
        scrollArea = new QScrollArea(centralwidget);
        scrollArea->setObjectName(QStringLiteral("scrollArea"));
        scrollArea->setGeometry(QRect(730, 10, 301, 951));
        scrollArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QStringLiteral("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 299, 949));
        groupBox_3 = new QGroupBox(scrollAreaWidgetContents);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        groupBox_3->setGeometry(QRect(10, 230, 281, 201));
        verticalLayoutWidget_3 = new QWidget(groupBox_3);
        verticalLayoutWidget_3->setObjectName(QStringLiteral("verticalLayoutWidget_3"));
        verticalLayoutWidget_3->setGeometry(QRect(10, 30, 261, 161));
        selectLayout1 = new QVBoxLayout(verticalLayoutWidget_3);
        selectLayout1->setObjectName(QStringLiteral("selectLayout1"));
        selectLayout1->setContentsMargins(0, 0, 0, 0);
        groupBox = new QGroupBox(scrollAreaWidgetContents);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        groupBox->setGeometry(QRect(10, 10, 281, 101));
        verticalLayoutWidget = new QWidget(groupBox);
        verticalLayoutWidget->setObjectName(QStringLiteral("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(9, 29, 261, 61));
        preprocessingLayout = new QVBoxLayout(verticalLayoutWidget);
        preprocessingLayout->setObjectName(QStringLiteral("preprocessingLayout"));
        preprocessingLayout->setContentsMargins(0, 0, 0, 0);
        groupBox_2 = new QGroupBox(scrollAreaWidgetContents);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        groupBox_2->setGeometry(QRect(10, 120, 281, 101));
        verticalLayoutWidget_2 = new QWidget(groupBox_2);
        verticalLayoutWidget_2->setObjectName(QStringLiteral("verticalLayoutWidget_2"));
        verticalLayoutWidget_2->setGeometry(QRect(10, 30, 261, 61));
        splitLayout = new QVBoxLayout(verticalLayoutWidget_2);
        splitLayout->setObjectName(QStringLiteral("splitLayout"));
        splitLayout->setContentsMargins(0, 0, 0, 0);
        groupBox_4 = new QGroupBox(scrollAreaWidgetContents);
        groupBox_4->setObjectName(QStringLiteral("groupBox_4"));
        groupBox_4->setGeometry(QRect(10, 440, 281, 101));
        verticalLayoutWidget_4 = new QWidget(groupBox_4);
        verticalLayoutWidget_4->setObjectName(QStringLiteral("verticalLayoutWidget_4"));
        verticalLayoutWidget_4->setGeometry(QRect(10, 30, 261, 61));
        mergeLayout = new QVBoxLayout(verticalLayoutWidget_4);
        mergeLayout->setObjectName(QStringLiteral("mergeLayout"));
        mergeLayout->setContentsMargins(0, 0, 0, 0);
        groupBox_5 = new QGroupBox(scrollAreaWidgetContents);
        groupBox_5->setObjectName(QStringLiteral("groupBox_5"));
        groupBox_5->setGeometry(QRect(10, 550, 281, 201));
        verticalLayoutWidget_5 = new QWidget(groupBox_5);
        verticalLayoutWidget_5->setObjectName(QStringLiteral("verticalLayoutWidget_5"));
        verticalLayoutWidget_5->setGeometry(QRect(10, 30, 261, 161));
        selectLayout2 = new QVBoxLayout(verticalLayoutWidget_5);
        selectLayout2->setObjectName(QStringLiteral("selectLayout2"));
        selectLayout2->setContentsMargins(0, 0, 0, 0);
        scrollArea->setWidget(scrollAreaWidgetContents);
        viewSelected2 = new QGraphicsView(centralwidget);
        viewSelected2->setObjectName(QStringLiteral("viewSelected2"));
        viewSelected2->setGeometry(QRect(0, 640, 361, 321));
        viewRectangles = new QGraphicsView(centralwidget);
        viewRectangles->setObjectName(QStringLiteral("viewRectangles"));
        viewRectangles->setGeometry(QRect(360, 640, 361, 321));
        SegmentWindow->setCentralWidget(centralwidget);
        statusBar = new QStatusBar(SegmentWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        SegmentWindow->setStatusBar(statusBar);

        retranslateUi(SegmentWindow);

        QMetaObject::connectSlotsByName(SegmentWindow);
    } // setupUi

    void retranslateUi(QMainWindow *SegmentWindow)
    {
        SegmentWindow->setWindowTitle(QApplication::translate("SegmentWindow", "MainWindow", 0));
        groupBox_3->setTitle(QApplication::translate("SegmentWindow", "Select 1 /2", 0));
        groupBox->setTitle(QApplication::translate("SegmentWindow", "Pre processing", 0));
        groupBox_2->setTitle(QApplication::translate("SegmentWindow", "Segment", 0));
        groupBox_4->setTitle(QApplication::translate("SegmentWindow", "Merge", 0));
        groupBox_5->setTitle(QApplication::translate("SegmentWindow", "Select 2 /2", 0));
    } // retranslateUi

};

namespace Ui {
    class SegmentWindow: public Ui_SegmentWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SEGMENTWINDOW_H
