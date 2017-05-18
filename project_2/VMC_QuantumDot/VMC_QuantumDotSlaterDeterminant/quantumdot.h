#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include "quantumstate.h"
#include <vector>
#include <string>
#include "particle.h"
#include <random>
#include "vec3.h"
#include <armadillo>

class QuantumDot
{
public:
    QuantumDot(double, int);
    int NumberOfParticles;
    double homega;
    double alpha;
    double beta;

    void getQuantumDotParticlesCoordinates();
    void getQuantumDotStates();
    void setVariationalParameters(double, double);
    void applyVMC(int);
    void applyVMCMPI(int);
    void setCoulombInterraction(int);
    void setJastrowFactor(int);
    double Alpha() { return alpha; }
    double Beta() { return beta; }
    void applySteepestDescent(int, int, double, double);

private:
    std::vector<Particle*> m_particles;
    std::vector<QuantumState> m_shells;
    vec3 m_QForceOld;
    vec3 m_QForceNew;

    int m_sm = -1;
    double m_sqrtomega;
    double m_alphaomega;
    arma::mat m_SD_down_inverse;
    arma::mat m_SD_up_inverse;
    arma::mat m_SD_down;
    arma::mat m_SD_up;
    double m_RSD; //Ratio of SD´s
    double m_RJ; //Ratio of Jartrow´s factors
    double m_gradientDotProductJastrowAllParticles;
    double m_LaplasianSD;
    double m_LaplasianJastrow;
    double m_DotProdGradientJastrowAndSD;
    const double m_s = 0.5;
    const double D = 0.5;     // difussion constant for importance sampling
    const double dt = 0.001; // dt for importance sampling
    int m_Coulomb;
    int m_Jastrow;
    double m_homega2;
    double m_alpha2;
    double m_alphaomega2;
    double m_KineticEnergy;
    double m_PotentialEnergy;
    double m_MeanRelativeDistance;

    typedef double (*polynomialArray) (double x); //reference to function type
    polynomialArray * polyRefArray = new polynomialArray[4]; //array with refs to Hermite polynomials
    typedef double (*FirstDerivativeArray) (double, double, double, double, double);
    FirstDerivativeArray * FirstDerivRefArrayX = new FirstDerivativeArray[4]; //array with refs to first derivatives X
    FirstDerivativeArray * FirstDerivRefArrayY = new FirstDerivativeArray[4]; //array with refs to first derivatives Y
    typedef double (*SecondDerivativeArray) (double, double, double, double, double);
    typedef double (*AlphaDerivativeArray) (double, double, double, double, double);
    AlphaDerivativeArray *AlphaDerivativeRefArray = new AlphaDerivativeArray[6]; // array with refs to first derivatives wrt alpha
    SecondDerivativeArray * SecondDerivRefArray = new SecondDerivativeArray[4]; //array with refs to second derivatives
    //Hermite polinomials
    static double H0(double x) { return 1.0; }
    static double H1(double x) { return 2.0*x; }
    static double H2(double x) { return 4.0*x*x - 2.0; }
    static double H3(double x) { return 8.0*x*x*x -12.0*x; }
    //gradient of the StateFunction X
    static double Fi0DerivativeX(double x, double y, double exponent, double aom, double sqrtom) { return -exponent*x*aom; }
    static double Fi1DerivativeX(double x, double y, double exponent, double aom, double sqrtom) { return -2.0*exponent*sqrtom*x*y*aom; }
    static double Fi2DerivativeX(double x, double y, double exponent, double aom, double sqrtom) { return 2.0*exponent*sqrtom*(1.0-x*x*aom); }
    static double Fi3DerivativeX(double x, double y, double exponent, double aom, double sqrtom) { return 1; }
    //gradient of the StateFunction Y
    static double Fi0DerivativeY(double x, double y, double exponent, double aom, double sqrtom) { return -exponent*y*aom; }
    static double Fi1DerivativeY(double x, double y, double exponent, double aom, double sqrtom) { return 2.0*exponent*sqrtom*(1.0-y*y*aom); }
    static double Fi2DerivativeY(double x, double y, double exponent, double aom, double sqrtom) { return -2.0*exponent*sqrtom*x*y*aom; }
    static double Fi3DerivativeY(double x, double y, double exponent, double aom, double sqrtom) { return 1; }
    //laplasian of the StateFunction
    static double Fi0Derivative2(double x, double y, double exponent, double aom, double sqrtom) { return exponent*aom*(aom*(x*x + y*y) - 2.0); }
    static double Fi1Derivative2(double x, double y, double exponent, double aom, double sqrtom) { return exponent*2.0*sqrtom*aom*y*(aom*(x*x + y*y)-4.0); }
    static double Fi2Derivative2(double x, double y, double exponent, double aom, double sqrtom) { return exponent*2.0*sqrtom*aom*x*(aom*(x*x + y*y)-4.0); }
    static double Fi3Derivative2(double x, double y, double exponent, double aom, double sqrtom) { return 1; }
    //Derivatives for steepest decent
    static double Fi0AlphaDerivative(double x, double y, double exponent, double sqrtom, double omega) { return -0.5*omega*(x*x+y*y)*exponent; }
    static double Fi1AlphaDerivative(double x, double y, double exponent, double sqrtom, double omega) { return -omega*(x*x+y*y)*exponent*sqrtom*y; }
    static double Fi2AlphaDerivative(double x, double y, double exponent, double sqrtom, double omega) { return -omega*(x*x+y*y)*exponent*sqrtom*x; }
    static double Fi3AlphaDerivative(double x, double y, double exponent, double sqrtom, double omega) { return -omega*(x*x+y*y)*exponent*(2*omega*y*y-1); }
    static double Fi4AlphaDerivative(double x, double y, double exponent, double sqrtom, double omega) { return -2*omega*omega*(x*x+y*y)*exponent*x*y; }
    static double Fi5AlphaDerivative(double x, double y, double exponent, double sqrtom, double omega) { return -omega*(x*x+y*y)*exponent*(2*omega*x*x-1); }

    void setUpStatesCartesian(int); //sets up nx/ny for cartesian basis
    void setUpCoordinatesCartesian(int); //sets up initial coordinates
    void setUpSlaterDeterminant(); //set initial SD
    polynomialArray getHetmitePolinomial(int); //input is a polin. order
    FirstDerivativeArray getGradientX(int);
    FirstDerivativeArray getGradientY(int);
    SecondDerivativeArray getLaplasian(int);
    AlphaDerivativeArray getDerivativeSDonAlpha(int); // input a state function number in SD

    double calculateSDRatio(size_t);
    double calculateGreenFunctionRatio(size_t);
    double calculateJastrowRatio(size_t);
    double calculateLocalEnergy();
    void calculateQuantumForce(size_t);
    void calculateQuantumForceNew(size_t);
    void updateInverseSlaterDeterminant(size_t);
    void updateSlaterDeterminant(size_t);
    void calculateDotProdGradientJastrowAndSD();
    void calculateLaplasianJastrow();
    void calculateLaplasianSD();
    void applyVMCSteepestDescent(int);
    void steepestDescentCalculateWFderivativeOnVarParameters();
    void resetSteepestDescentHelpVars();
    void writeVectorToBinaryFile(std::string, std::vector<double>&);

    //SteepestDescent variables
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

};

#endif // QUANTUMDOT_H
