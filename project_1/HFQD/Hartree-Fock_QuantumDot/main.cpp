#include <iostream>
#include "quantumdot.h"
#include "Coulomb_Functions.hpp"

using namespace std;

int main()
{

    QuantumDot qdot(4);
    qdot.getQuantumDotStates();
    Coulomb_HO(double &hw, int &ni, int &mi, int &nj, int &mj, int &nk, int &mk, int &nl, int &ml)
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

