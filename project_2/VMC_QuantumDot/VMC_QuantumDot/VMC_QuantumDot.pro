TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    quantumstate.cpp \
    quantumdot.cpp \
    vec3.cpp \
    particle.cpp \
    quantumforce.cpp

HEADERS += \
    quantumstate.h \
    quantumdot.h \
    vec3.h \
    particle.h \
    random.h \
    quantumforce.h
