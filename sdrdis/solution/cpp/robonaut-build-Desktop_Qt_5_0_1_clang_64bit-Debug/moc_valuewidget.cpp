/****************************************************************************
** Meta object code from reading C++ file 'valuewidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.0.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../robonaut/common/valuewidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'valuewidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.0.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_ValueWidget_t {
    QByteArrayData data[9];
    char stringdata[83];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    offsetof(qt_meta_stringdata_ValueWidget_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData) \
    )
static const qt_meta_stringdata_ValueWidget_t qt_meta_stringdata_ValueWidget = {
    {
QT_MOC_LITERAL(0, 0, 11),
QT_MOC_LITERAL(1, 12, 22),
QT_MOC_LITERAL(2, 35, 0),
QT_MOC_LITERAL(3, 36, 5),
QT_MOC_LITERAL(4, 42, 19),
QT_MOC_LITERAL(5, 62, 4),
QT_MOC_LITERAL(6, 67, 3),
QT_MOC_LITERAL(7, 71, 3),
QT_MOC_LITERAL(8, 75, 6)
    },
    "ValueWidget\0on_slider_valueChanged\0\0"
    "value\0on_text_textChanged\0arg1\0min\0"
    "max\0nbStep\0"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ValueWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       4,   30, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   24,    2, 0x08,
       4,    1,   27,    2, 0x08,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::QString,    5,

 // properties: name, type, flags
       3, QMetaType::Float, 0x00095103,
       6, QMetaType::Float, 0x00095103,
       7, QMetaType::Float, 0x00095103,
       8, QMetaType::Int, 0x00095103,

       0        // eod
};

void ValueWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ValueWidget *_t = static_cast<ValueWidget *>(_o);
        switch (_id) {
        case 0: _t->on_slider_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->on_text_textChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject ValueWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ValueWidget.data,
      qt_meta_data_ValueWidget,  qt_static_metacall, 0, 0}
};


const QMetaObject *ValueWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ValueWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ValueWidget.stringdata))
        return static_cast<void*>(const_cast< ValueWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int ValueWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 2)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 2;
    }
#ifndef QT_NO_PROPERTIES
      else if (_c == QMetaObject::ReadProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< float*>(_v) = value(); break;
        case 1: *reinterpret_cast< float*>(_v) = min(); break;
        case 2: *reinterpret_cast< float*>(_v) = max(); break;
        case 3: *reinterpret_cast< int*>(_v) = nbStep(); break;
        }
        _id -= 4;
    } else if (_c == QMetaObject::WriteProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: setValue(*reinterpret_cast< float*>(_v)); break;
        case 1: setMin(*reinterpret_cast< float*>(_v)); break;
        case 2: setMax(*reinterpret_cast< float*>(_v)); break;
        case 3: setNbStep(*reinterpret_cast< int*>(_v)); break;
        }
        _id -= 4;
    } else if (_c == QMetaObject::ResetProperty) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 4;
    } else if (_c == QMetaObject::RegisterPropertyMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}
QT_END_MOC_NAMESPACE
