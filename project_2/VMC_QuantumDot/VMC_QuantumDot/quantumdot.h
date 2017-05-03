#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include "quantumstate.h"
#include <vector>
#include "particle.h"
#include <random>
#include "quantumforce.h"

class QuantumDot
{
public:
    QuantumDot(double, int);
    int NumberOfParticles;
    double homega;
    double alpha;
    double beta;
    void setVariationalParameters(double, double);
    void getQuantumDotStates();
    void getQuantumDotStatesNumber();
    void getQuantumDotParticlesCoordinates();
    void applyVMC(int);
    void applyVMCstandard(int);
    double Alpha() { return alpha; };
    double Beta() { return beta; };
    void setCoulombInterraction(int);
    void setJastrowFactor(int);
    void setSpinParameter(int);

private:
    std::vector<Particle*> m_particles;
    std::vector<QuantumState> m_shells;
    std::vector<double> m_localEnergy;
    void initialize(int);
    int m_sm = -1;
    const double m_s = 0.5;
    const double D = 0.5;     // difussion constant for importance sampling
    const double dt = 0.01; // dt for importance sampling
    int m_Coulomb;
    int m_Jastrow;
    double m_homega2;
    double m_alpha2;
    double m_alphaomega;
    double m_alphaomega2;
    double m_a;
    void setUpStatesCartesian(int);
    double calculateLocalEnergy();
    double calculateLocalEnergyWithoutJastrow();
    double calculateGreenFunctionRatio(size_t);
    double calculateJastrowRatio(size_t);
    void writeLocalEnergyToFile(double, string);
    void writeVectorToFile(string);


    //double calculateKineticEnergyNumerical();
    //double m_RelativeDistanceOld;
    //double m_RelativeDistanceNew;
    //double computeRelativeDistance();


};

#endif // QUANTUMDOT_H
