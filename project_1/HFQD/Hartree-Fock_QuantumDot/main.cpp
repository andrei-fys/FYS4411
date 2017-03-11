#include <iostream>
#include "quantumdot.h"
//#include "Coulomb_Functions.hpp"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int NumberOfShells = 3;
    int NumberOfElectrons = 6;
    int HOStrenth = 1;


    // If a first argument is provided, it is the number of shells
    if(numberOfArguments > 1) NumberOfShells = atoi(argumentList[1]);
    // If a second argument is provided, it is the number of electrons
    if(numberOfArguments > 2) NumberOfElectrons = atoi(argumentList[2]);
    // If a third argument is provided, it is the HO strenth /omega
    if(numberOfArguments > 3) HOStrenth = atoi(argumentList[3]);

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

