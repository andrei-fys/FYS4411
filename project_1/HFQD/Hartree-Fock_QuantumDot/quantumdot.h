#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include <vector>
#include <iostream>
#include "quantumstate.h"

class QuantumDot
{
public:
    QuantumDot(int);
    int EnergyCutOff;
    void getQuantumDotStates();
private:
    std::vector<QuantumState> m_shells;
    int m_sm = -1;
    const double m_s = 0.5;
    void setUpStates(int);
};

#endif // QUANTUMDOT_H
