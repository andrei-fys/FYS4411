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

void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << quantum_state.nx() << endl;
        cout << quantum_state.ny() << endl;
        cout << quantum_state.sm() << endl;
        cout << quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

