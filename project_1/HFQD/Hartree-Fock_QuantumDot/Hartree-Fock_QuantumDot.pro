TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    quantumdot.cpp \
    quantumstate.cpp \
    Coulomb_Functions.cpp

HEADERS += \
    quantumdot.h \
    quantumstate.h \
    Coulomb_Functions.hpp

