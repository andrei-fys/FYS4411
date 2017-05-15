TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -O3

INCLUDEPATH += /usr/local/include
INCLUDEPATH += ~/armadillo_install/include
LIBS += -L/usr/local/lib
LIBS += -L~/armadillo_install/lib
LIBS += -llapack -lblas -larmadillo
##Uncomment for MPI support
##MPI Settings
#QMAKE_CXX = mpicxx
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = mpicc
#
#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK



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
