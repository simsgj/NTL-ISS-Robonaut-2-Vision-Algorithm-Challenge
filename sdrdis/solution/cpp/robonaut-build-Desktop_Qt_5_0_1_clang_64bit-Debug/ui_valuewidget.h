/********************************************************************************
** Form generated from reading UI file 'valuewidget.ui'
**
** Created by: Qt User Interface Compiler version 5.0.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VALUEWIDGET_H
#define UI_VALUEWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QSlider>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ValueWidget
{
public:
    QSlider *slider;
    QLineEdit *text;
    QLabel *label;

    void setupUi(QWidget *ValueWidget)
    {
        if (ValueWidget->objectName().isEmpty())
            ValueWidget->setObjectName(QStringLiteral("ValueWidget"));
        ValueWidget->resize(411, 21);
        slider = new QSlider(ValueWidget);
        slider->setObjectName(QStringLiteral("slider"));
        slider->setGeometry(QRect(110, 0, 251, 22));
        slider->setMinimum(0);
        slider->setOrientation(Qt::Horizontal);
        slider->setInvertedAppearance(false);
        slider->setInvertedControls(false);
        text = new QLineEdit(ValueWidget);
        text->setObjectName(QStringLiteral("text"));
        text->setGeometry(QRect(370, 0, 41, 21));
        label = new QLabel(ValueWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(0, 0, 101, 21));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        retranslateUi(ValueWidget);

        QMetaObject::connectSlotsByName(ValueWidget);
    } // setupUi

    void retranslateUi(QWidget *ValueWidget)
    {
        ValueWidget->setWindowTitle(QApplication::translate("ValueWidget", "ValueWidget", 0));
        label->setText(QApplication::translate("ValueWidget", "TextLabel", 0));
    } // retranslateUi

};

namespace Ui {
    class ValueWidget: public Ui_ValueWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VALUEWIDGET_H
