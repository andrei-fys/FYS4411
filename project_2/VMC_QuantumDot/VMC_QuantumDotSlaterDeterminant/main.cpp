#include <iostream>
#include "quantumdot.h"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int NumberOfElectrons = 12;
    double HOStrenth = 1.0;
    int MonteCarloSamples = (int) 1e7;
    double alpha = 0.624204;   //if steepest descent is used, values are initial values for var. params
    double beta = 2.01695;

    int MonteCarloSamplesVariational = 1000000;
    int MaxSteepestDescentIterations = 200;
    double SteepestDescentStep = 0.001;
    double tolerance = 10e-5;

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

    qdot.applySteepestDescent(MonteCarloSamplesVariational,
                              MaxSteepestDescentIterations,
                              SteepestDescentStep,
                              tolerance);

    //qdot.applyVMC(MonteCarloSamples);
    //qdot.getQuantumDotStates();
    //qdot.getQuantumDotParticlesCoordinates();

}
