#include "quantumstate.h"

QuantumState::QuantumState()
{

}

void QuantumState::set(int nx, int ny, int sm, double s){
    m_nx = nx;
    m_ny = ny;
    m_sm = sm;
    m_s = s;
}

void QuantumState::flipSpin(){
    m_sm = m_sm*(-1);
}
