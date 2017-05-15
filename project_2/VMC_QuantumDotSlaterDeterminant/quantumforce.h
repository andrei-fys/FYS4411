#ifndef QUANTUMFORCE_H
#define QUANTUMFORCE_H


class QuantumForce
{
public:
    QuantumForce(double, double);
    double evaluate(double, double, double, double);
private:
    double m_omega;
    double m_alpha;
};

#endif // QUANTUMFORCE_H
