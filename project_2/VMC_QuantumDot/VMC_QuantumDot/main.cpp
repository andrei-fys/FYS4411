#include <iostream>
#include "quantumdot.h"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int NumberOfElectrons = 2;
    double HOStrenth = 1;
    //int MonteCarloSamples = 1048576; //2^20 = 1048576;
    int MonteCarloSamples = 8388608; //2¨^23
    double alpha = 1.0;
    double beta = 0.0;


    // If a first argument is provided, it is the number of electrons
    if(numberOfArguments > 1) NumberOfElectrons = atoi(argumentList[1]);
    // If a second argument is provided, it is the HO strenth /omega
    if(numberOfArguments > 2) HOStrenth = atof(argumentList[2]);
    // If a third argument is provided, it is the number of Monte Carlo samples
    if(numberOfArguments > 2) MonteCarloSamples = atof(argumentList[3]);


    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setCoulombInterraction(1); // 0 - to turn off Coulomb interraction
    qdot.setVariationalParameters(alpha, beta);
    qdot.applyVMC(MonteCarloSamples);
    cout << "=================================" << endl;
    qdot.applyVMCstandard(MonteCarloSamples);

    //qdot.getQuantumDotParticlesCoordinates();


}
