#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include <vector>
#include <iostream>
#include "quantumstate.h"
#include "Coulomb_Functions.hpp"
#include <armadillo>

class QuantumDot{
public:
    QuantumDot(int, double, int);
    int EnergyCutOff;
    int NumberOfParticles;
    double homega;
    void getQuantumDotStates();
    //arma::mat computeDensityMatrix(arma::mat);

    void diagonalizeHFMatrix();
    void getQuantumDotStatesNumber();
    void applyHartreeFockMethod();
private:
    std::vector<QuantumState> m_shells;
    int m_sm = -1;
    const double m_s = 0.5;
    arma::mat m_C;              //Ciefficient matrix
    //arma::mat m_DensityMatrix;
    arma::mat m_HF;             //Hartree-Fock matrix
    arma::vec eigval_previous;
    arma::vec eigval;
    arma::mat eigvec;
    void setUpStatesCartesian(int);
    void setUpStatesPolar(int, double h_omega, int);
    void setUpStatesPolarSorted(int, double h_omega, int);
    void setCoefficientMatrix(arma::mat);
    void computeHFmatrix(arma::mat);
    arma::mat computeDensityMatrix();
    double CalculateNonIntEnergy(int n, int m);
    double computeHartreeFockEnergyDifference();
    void computeHartreeFockEnergy();

};

#endif // QUANTUMDOT_H
