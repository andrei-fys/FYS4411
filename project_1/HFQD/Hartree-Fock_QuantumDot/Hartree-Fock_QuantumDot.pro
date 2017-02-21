TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    quantumdot.cpp \
    quantumstate.cpp \
    Coulomb_Functions.cpp \
    QD_Coulomb.cpp

HEADERS += \
    quantumdot.h \
    quantumstate.h \
    Coulomb_Functions.hpp

