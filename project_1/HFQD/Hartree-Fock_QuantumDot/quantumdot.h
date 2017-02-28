#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include <vector>
#include <iostream>
#include "quantumstate.h"
#include "Coulomb_Functions.hpp"

class QuantumDot
{
public:
    QuantumDot(int);
    int EnergyCutOff;
    void getQuantumDotStates();
    double computeCoulombInteractionPolar();
    void getQuantumDotStatesNumber();
private:
    std::vector<QuantumState> m_shells;
    int m_sm = -1;
    const double m_s = 0.5;
    void setUpStatesCartesian(int);
    void setUpStatesPolar(int);
};

#endif // QUANTUMDOT_H
