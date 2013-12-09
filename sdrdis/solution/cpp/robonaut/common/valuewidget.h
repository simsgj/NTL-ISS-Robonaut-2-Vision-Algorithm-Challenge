#ifndef VALUEWIDGET_H
#define VALUEWIDGET_H

#include <QWidget>
#include <QLabel>
#include <QDebug>


#ifdef _DESIGNER_EXPORT
#include <QDesignerExportWidget>
#else
#define QDESIGNER_WIDGET_EXPORT
#endif


namespace Ui {
class ValueWidget;
}

class QDESIGNER_WIDGET_EXPORT ValueWidget : public QWidget
{
    Q_OBJECT

public:

    Q_PROPERTY(float value READ value WRITE setValue)
    Q_PROPERTY(float min READ min WRITE setMin)
    Q_PROPERTY(float max READ max WRITE setMax)
    Q_PROPERTY(int nbStep READ nbStep WRITE setNbStep)



    explicit ValueWidget(QWidget *parent = 0);
    ~ValueWidget();
    void refresh();

    void setValue(float value);
    float value();

    void setMin(float min);
    float min();

    void setMax(float max);
    float max();

    void setNbStep(int nbStep);
    int nbStep();

private slots:
    void on_slider_valueChanged(int value);

    void on_text_textChanged(const QString &arg1);

private:
    Ui::ValueWidget *ui;
    QString labelName;
    float _value;
    float _min;
    float _max;
    int _nbStep;
};

#endif // VALUEWIDGET_H
