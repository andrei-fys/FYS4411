#include "quantumforce.h"

QuantumForce::QuantumForce(double alpha, double omega){
    m_omega = omega;
    m_alpha = alpha;
}


double QuantumForce::evaluate(double x1, double y1, double x2, double y2){
    return -2.0*m_alpha*m_omega*(x1 + y1 + x2 + y2);
}
