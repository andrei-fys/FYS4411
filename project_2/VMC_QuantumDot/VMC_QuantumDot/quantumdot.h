#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include "quantumstate.h"
#include <vector>
#include "particle.h"
#include <random>
#include "vec3.h"
//#include "quantumforce.h"

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
    double Alpha() { return alpha; }
    double Beta() { return beta; }
    void setCoulombInterraction(int);
    void setJastrowFactor(int);
    void setSpinParameter(int);
    void applySteepestDescent(int, int, double, double);

private:
    std::vector<Particle*> m_particles;
    std::vector<QuantumState> m_shells;
    std::vector<double> m_localEnergy;
    vec3 m_QForceOld;
    vec3 m_QForceNew;
    void initialize(int);
    int m_sm = -1;
    const double m_s = 0.5;
    const double D = 0.5;     // difussion constant for importance sampling
    const double dt = 0.001; // dt for importance sampling
    int m_Coulomb;
    int m_Jastrow;
    double m_homega2;
    double m_alpha2;
    double m_alphaomega;
    double m_alphaomega2;
    double m_a;
    //double m_Energy;
    //SteepestDescent variables
    bool m_InsideSteepestDescent;
    double m_ExpectationLocalEnergyDerivativeAlphaSecondTerm;
    double m_ExpectationLocalEnergyDerivativeBetaSecondTerm;
    double m_ExpectationLocalEnergyDerivativeAlphaFirstTerm;
    double m_ExpectationLocalEnergyDerivativeBetaFirstTerm;
    double m_LocalEnergyWFDerivativeAlpha;
    double m_LocalEnergyWFDerivativeBeta;
    double m_MeanLocalEnergyWFDerivativeAlpha;
    double m_MeanLocalEnergyWFDerivativeBeta;
    double m_ExpectationLocalEnergyDerivativeAlpha;
    double m_ExpectationLocalEnergyDerivativeBeta;

    void setUpStatesCartesian(int);
    double calculateLocalEnergy();
    double calculateLocalEnergyWithoutJastrow();
    double calculateGreenFunctionRatio(size_t);
    double calculateJastrowRatio(size_t);
    void writeLocalEnergyToFile(double, string);
    void writeVectorToFile(string);
    void resetSteepestDescentHelpVars();
    void calculateQuantumForce(size_t);
    void calculateQuantumForceNew(size_t);



    //double calculateKineticEnergyNumerical();
    //double m_RelativeDistanceOld;
    //double m_RelativeDistanceNew;
    //double computeRelativeDistance();


};

#endif // QUANTUMDOT_H
