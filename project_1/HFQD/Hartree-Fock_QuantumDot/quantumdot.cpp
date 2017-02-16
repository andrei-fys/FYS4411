#include "quantumdot.h"
#include <iostream>
#include <vector>

using namespace std;

QuantumDot::QuantumDot(int EnergyCutOff){
    setUpStates(EnergyCutOff);

}


void QuantumDot::setUpStates(int EnergyCutOff) {
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

void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << quantum_state.nx() << endl;
        cout << quantum_state.ny() << endl;
        cout << quantum_state.sm() << endl;
        cout << quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

