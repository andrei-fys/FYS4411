#include "quantumdot.h"
#include <iostream>
#include <vector>

using namespace std;

QuantumDot::QuantumDot(int EnergyCutOff){
    //setUpStatesCartesian(EnergyCutOff);
    setUpStatesPolar(EnergyCutOff);

}


void QuantumDot::setUpStatesCartesian(int EnergyCutOff) {
    //QuantumState * m_q_state = new QuantumState;
    QuantumState m_q_state;
    int nx, ny;
    for (int i = 1; i < EnergyCutOff + 1; i++){
        for (int j = 0; j < i; j++){
            nx = j;
            ny = i - nx -1;
            //QuantumState * m_q_state = new QuantumState;
            m_q_state.set(nx, ny, m_sm, m_s);
            m_shells.push_back(m_q_state);
            m_q_state.flipSpin();
            m_shells.push_back(m_q_state);
         }
    }
}

void QuantumDot::setUpStatesPolar(int EnergyCutOff) {
    //QuantumState * m_q_state = new QuantumState;
    QuantumState m_q_state;
    vector<int> oddShells;
    int n, m;
    //loop to find all odd shell numbers including energyCutOff shell
    for (int i = 1; i <= EnergyCutOff; i++){
        if (i % 2 != 0){
            oddShells.push_back(i);
        }
    }
    for (int index = 0; index < oddShells.size(); ++index)
    {
       //cout << index << "  " << oddShells[index] << endl;
       //positive m
       for (int j=0; j<=EnergyCutOff - oddShells[index]; j++){
           //cout << "[ " << index << " , " << j << " ]" << endl;
           n = index;
           m = j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
           m_q_state.flipSpin();
           m_shells.push_back(m_q_state);
       }
       //negative m
       for (int j=1; j<=EnergyCutOff - oddShells[index]; j++){
           //cout << "[ " << index << " , " << j << " ]" << endl;
           n = index;
           m = -1*j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
           m_q_state.flipSpin();
           m_shells.push_back(m_q_state);
       }
    }

}

double QuantumDot::computeCoulombInteractionPolar(){
    for(QuantumState quantum_state : m_shells){
        int p_n = quantum_state.n();
        int p_m = quantum_state.m();
        int p_sm = quantum_state.sm();
        int p_s = quantum_state.s();
        for(QuantumState quantum_state : m_shells){
            int q_n = quantum_state.n();
            int q_m = quantum_state.m();
            int q_sm = quantum_state.sm();
            int q_s = quantum_state.s();
            for(QuantumState quantum_state : m_shells){
                int r_n = quantum_state.n();
                int r_m = quantum_state.m();
                int r_sm = quantum_state.sm();
                int r_s = quantum_state.s();
                for(QuantumState quantum_state : m_shells){
                    int s_n = quantum_state.n();
                    int s_m = quantum_state.m();
                    int s_sm = quantum_state.sm();
                    int s_s = quantum_state.s();
                    if (p_m + q_m == r_m + s_m){
                        if (p_sm + q_sm == r_sm + s_sm){
                            if (p_sm == r_sm && q_sm == s_sm){
                                //double TBME = Coulomb_HO(1.0, p_n, p_sm, q_n, q_sm, r_n, r_sm, s_n, s_sm);
                                //cout << setprecision(12);
                                //cout << "< " << p_n << "," << p_sm << " ; " << q_n << "," << q_sm << " || V || " << r_n << "," << r_sm << " ; " << s_n << "," << s_sm << " > = " << TBME << std::endl;
                                cout << "< " << p_n << "," << p_sm << " ; " << q_n << "," << q_sm << " || V || " << r_n << "," << r_sm << " ; " << s_n << "," << s_sm << " > = " << endl;
                            } else {
                                continue;
                                //cout << "third check 0 " << endl;
                            }
                        } else {
                            continue;
                            //cout << "second check 0 " << endl;
                        }
                    } else {
                        continue;
                        //cout << "first check 0 " << endl;
                    }
                }
            }
        }

    }
return 12.3;
}

void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << quantum_state.n() << endl;
        cout << quantum_state.m() << endl;
        cout << quantum_state.sm() << endl;
        cout << quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

void QuantumDot::getQuantumDotStatesNumber(){
    cout << m_shells.size() << endl;
}

