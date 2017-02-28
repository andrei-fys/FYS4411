#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H


class QuantumState
{
private:
    int m_n;     // n
    int m_m;     // m
    int m_sm;     // spin projection (1 or -1)
    double m_s;   // spin

public:
    QuantumState();
    double s() { return m_s; }
    int n() { return m_n; }
    int m() { return m_m; }
    int sm() { return m_sm; }
    void set(int n, int m, int sm, double s);
    void flipSpin();
};

#endif // QUANTUMSTATE_H
