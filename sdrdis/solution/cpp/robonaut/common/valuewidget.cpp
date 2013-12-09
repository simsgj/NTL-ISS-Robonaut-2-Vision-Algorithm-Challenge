#include "valuewidget.h"
#include "ui_valuewidget.h"

ValueWidget::ValueWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ValueWidget)
{
    ui->setupUi(this);
    this->_min = -1;
    this->_max = 1;
    this->_value = 0;
    this->setNbStep(100);
    this->labelName = "Label: ";
    this->refresh();
}

ValueWidget::~ValueWidget()
{
    delete ui;
}

void ValueWidget::on_slider_valueChanged(int value)
{
    this->setValue((value * 1.0 / this->_nbStep) * (this->_max - this->_min) + this->_min);
}

void ValueWidget::refresh() {
    bool *ok = new bool(true);
    float value = ui->text->text().toFloat(ok);
    if ((*ok) == false || (value - this->_value) > 1e-4 || (value - this->_value) < -1e-4) {
        ui->text->setText(QString().setNum(this->_value));
    }
    delete ok;
    ui->slider->setValue(((this->_value - this->_min) / (this->_max - this->_min)) * this->_nbStep);
    ui->label->setText(this->labelName);
}

void ValueWidget::on_text_textChanged(const QString &arg1)
{
    bool *ok = new bool(true);
    float value = arg1.toFloat(ok);
    if ((*ok) == true) {
        this->setValue(value);
    }
    delete ok;
}

void ValueWidget::setValue(float value) {
    if (this->_value < this->_min) {
        this->_value = this->_min;
    }
    if (this->_value > this->_max) {
        this->_value = this->_max;
    }
    this->_value = value;
    this->refresh();
}

float ValueWidget::value() {
    return this->_value;
}

void ValueWidget::setMin(float min) {
    this->_min = min;
    this->refresh();
}

float ValueWidget::min() {
    return this->_min;
}

void ValueWidget::setMax(float max) {
    this->_max = max;
    this->refresh();
}

float ValueWidget::max() {
    return this->_max;
}

void ValueWidget::setNbStep(int nbStep) {
    this->_nbStep = nbStep;
    this->ui->slider->setMaximum(nbStep);
    this->refresh();
}

int ValueWidget::nbStep() {
    return this->_nbStep;
}
