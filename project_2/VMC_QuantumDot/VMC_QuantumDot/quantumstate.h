#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H


class QuantumState
{
private:
    int m_nx;     // nx
    int m_ny;     // ny
    int m_sm;     // spin projection (1 or -1)
    double m_s;   // spin

public:
    QuantumState();
    double s() { return m_s; }
    int nx() { return m_nx; }
    int ny() { return m_ny; }
    int sm() { return m_sm; }
    void set(int nx, int ny, int sm, double s);
    void flipSpin();
};

#endif // QUANTUMSTATE_H
