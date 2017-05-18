#include <iostream>
#include "quantumdot.h"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int NumberOfElectrons = 6;
    double HOStrenth = 1.0;
    int MonteCarloSamples = (int) 1e6;//1000000;
    double alpha = 1;
    double beta = 0.4;

    int MonteCarloSamplesVariational = 100000;
    int MaxSteepestDescentIterations = 20;
    double SteepestDescentStep = 0.001;
    double tolerance = 0.7;

    // If a first argument is provided, it is the number of electrons
    if(numberOfArguments > 1) NumberOfElectrons = atoi(argumentList[1]);
    // If a second argument is provided, it is the HO strenth /omega
    if(numberOfArguments > 2) HOStrenth = atof(argumentList[2]);
    // If a third argument is provided, it is the number of Monte Carlo samples
    if(numberOfArguments > 2) MonteCarloSamples = atof(argumentList[3]);

    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setVariationalParameters(alpha, beta);
    qdot.setCoulombInterraction(1);
    qdot.setJastrowFactor(1);

    /*qdot.applySteepestDescent(MonteCarloSamplesVariational,
                              MaxSteepestDescentIterations,
                              SteepestDescentStep,
                              tolerance);
    */
    qdot.applyVMC(MonteCarloSamples);
    //qdot.getQuantumDotStates();
    //qdot.getQuantumDotParticlesCoordinates();

}
