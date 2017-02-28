#include "quantumstate.h"

QuantumState::QuantumState()
{

}

void QuantumState::set(int n, int m, int sm, double s){
    m_n = n;
    m_m = m;
    m_sm = sm;
    m_s = s;
}

void QuantumState::flipSpin(){
    m_sm = m_sm*(-1);
}
