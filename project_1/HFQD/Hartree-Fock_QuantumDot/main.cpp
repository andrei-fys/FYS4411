#include <iostream>
#include "quantumdot.h"
#include "Coulomb_Functions.hpp"

using namespace std;

int main()
{

    int NumberOfShells = 5;
    int HOStrenth = 1;
    int NumberOfElectrons = 6;
    QuantumDot qdot(NumberOfShells, HOStrenth, NumberOfElectrons);
    //qdot.getQuantumDotStates();


    qdot.applyHartreeFockMethod();

}



/*class QuatumState;

class MatrixElements {
    double computeOnebodyElement(i,j);
    double computeTwobodyElementAntiSymmetric(i,j,k,l);
};

class HartreeFock {
private:
    mat FockMatrix;
    mat Coefficients; // C
    mat DensityMatrix; // rho
    mat HartreeFockEnergies;
    mat epsilon;

public:
    mat computeFockMatrix();
    mat computeDensityMatrix();
    void diagonalizeFockMatrix();
    double computeHartreeoFockEnergy();

    void selfConsisintentFIeldIterations();
};
*/

