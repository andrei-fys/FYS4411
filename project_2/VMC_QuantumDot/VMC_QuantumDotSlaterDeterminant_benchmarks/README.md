By default MPI version is compiled. You need to have MPI installed.
To compile with QTCreator you need to have following in pro-file:


```bash
#Uncomment for MPI support
#MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile)
-DMPICH_IGNORE_CXX_SEEK
```
Additionally you need armadillo library. In pro-file should be following lines:

```bash
LIBS += -llapack -lblas -larmadillo
```
Run compiled version simply with:
```bash
mpirun -n 4 ./VMC_QuantumDotSD
```
It will give following output and writes file [output for 6 electrons](LocalEnergy_1.000000_6_el_0_results):
```bash
alpha 0.701856
beta 1.70348
Energy 20.4851
Variance 1.32016e-06
Kinetic 3.34429
Potential 17.1408
Mean RelDist 6.58931
Accept % :99.9873
Total number of MC samples 4000000
```
