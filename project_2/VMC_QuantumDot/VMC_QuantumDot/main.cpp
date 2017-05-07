#include <iostream>
#include "quantumdot.h"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int g_argc;
    char **g_argv;
    g_argc = numberOfArguments;
    g_argv = argumentList;

    int NumberOfElectrons = 2;
    double HOStrenth = 0.5;
    int MonteCarloSamples = 100000000;
    //int MonteCarloSamples = 8388608; //2^23
    int MonteCarloSamplesVariational = 1000000;
    double tolerance = 10e-7;
    int MaxSteepestDescentIterations = 20;
    double SteepestDescentStep = 0.05;
    double alpha = 0.989428;
    double beta = 0.309102;


    // If a first argument is provided, it is the number of electrons
    if(numberOfArguments > 1) NumberOfElectrons = atoi(argumentList[1]);
    // If a second argument is provided, it is the HO strenth /omega
    if(numberOfArguments > 2) HOStrenth = atof(argumentList[2]);
    // If a third argument is provided, it is the number of Monte Carlo samples
    if(numberOfArguments > 2) MonteCarloSamples = atof(argumentList[3]);

/*
    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setCoulombInterraction(1); // 0 - to turn off Coulomb interraction
    qdot.setJastrowFactor(1);       // 0 - to turn off correlations
    qdot.setSpinParameter(1);       // 1 - antiparallel, 1/3 - parallel
    qdot.setVariationalParameters(alpha, beta);
    qdot.applyVMC(MonteCarloSamples);
    cout << "=================================" << endl;
    //qdot.applyVMCstandard(MonteCarloSamples); //works without Jastrow factor, just with Coulomb on/off
    //qdot.getQuantumDotParticlesCoordinates();

*/

    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setCoulombInterraction(1); // 0 - to turn off Coulomb interraction
    qdot.setJastrowFactor(1);       // 0 - to turn off correlations
    qdot.setSpinParameter(1);       // 1 - antiparallel, 1/3 - parallel
    qdot.setVariationalParameters(alpha, beta);

    qdot.applySteepestDescent(MonteCarloSamplesVariational,
                              MaxSteepestDescentIterations,
                              SteepestDescentStep,
                              tolerance);

    qdot.applyVMC(MonteCarloSamples);
    cout << "=================================" << endl;


/*
    //MPI test
    alpha = 0.986561;
    beta = 0.415157;
    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setMPIenv(g_argc, g_argv);
    qdot.setCoulombInterraction(1); // 0 - to turn off Coulomb interraction
    qdot.setJastrowFactor(1);       // 0 - to turn off correlations
    qdot.setSpinParameter(1);       // 1 - antiparallel, 1/3 - parallel
    qdot.setVariationalParameters(alpha, beta);
    qdot.applyVMC(MonteCarloSamples);
/*

/*
    alpha 0.986561
    beta 0.415157
    Energy 3.00833
    Variance 8.63128e-11
    Accept 98.4854
*/

/*  // Without importance sampling
    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setCoulombInterraction(1); // 0 - to turn off Coulomb interraction
    qdot.setJastrowFactor(1);       // 0 - to turn off correlations
    qdot.setSpinParameter(1);       // 1 - antiparallel, 1/3 - parallel
    qdot.setVariationalParameters(alpha, beta);
    qdot.applyVMCstandard(MonteCarloSamples);
*/
}
