#include <iostream>
#include "quantumdot.h"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    clock_t start, finish;
    start = clock();
    int numprocs;
    int NumberOfElectrons = 6;
    double HOStrenth = 1;
    int MonteCarloSamples = (long int) 1e6;
    double alpha = 0.701856;   //if steepest descent is used, values are initial values for var. params
    double beta = 1.70348;

    int MonteCarloSamplesVariational = 1000000;
    int MaxSteepestDescentIterations = 200;
    int RunSteepestDescent = 0;
    double SteepestDescentStep = 0.1;
    double tolerance = 10e-8;

    // If a first argument is provided, it is the number of electrons
    if(numberOfArguments > 1) NumberOfElectrons = atoi(argumentList[1]);
    // If a second argument is provided, it is the HO strenth /omega
    if(numberOfArguments > 2) HOStrenth = atof(argumentList[2]);
    // If a third argument is provided, it is the number of Monte Carlo samples
    if(numberOfArguments > 3) MonteCarloSamples = atof(argumentList[3]);
    // If a fourth argument is provided, it is the first variational parameter - alpha
    if(numberOfArguments > 4) alpha = atof(argumentList[4]);
    // If a fivth argument is provided, it is the second variational parameter - beta
    if(numberOfArguments > 5) beta = atof(argumentList[5]);
    // If a sixth argument is provided, it is the Steepest Descent trigger, 0/1 (1 - to run)
    if(numberOfArguments > 6) RunSteepestDescent = atof(argumentList[6]);
    // If a seventh argument is provided, it is the number of Monte Carlo samples for Steepest Descent
    if(numberOfArguments > 7) MonteCarloSamplesVariational = atof(argumentList[7]);
    // If a eighth argument is provided, it is the Steepest Descent initial step
    if(numberOfArguments > 8) SteepestDescentStep = atof(argumentList[8]);

    QuantumDot qdot(HOStrenth, NumberOfElectrons);
    qdot.setVariationalParameters(alpha, beta);
    qdot.setCoulombInterraction(1);
    qdot.setJastrowFactor(1);

    if (RunSteepestDescent == 1) {
        qdot.applySteepestDescent(MonteCarloSamplesVariational,
                                  MaxSteepestDescentIterations,
                                  SteepestDescentStep,
                                  tolerance);
    }


    //qdot.applyVMC(MonteCarloSamples);


    qdot.applyVMCMPI(MonteCarloSamples, numprocs);
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << "Time total = " << timeused  << endl;
}
