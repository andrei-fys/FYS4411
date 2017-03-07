#include "quantumdot.h"
#include <iostream>
//#include <vector>
//#include <armadillo>

using namespace std;

bool myComparison(const pair<int,double> &a,const pair<int,double> &b)
{
       return a.second<b.second;
}

QuantumDot::QuantumDot(int EnergyCutOff, double h_omega, int ParticlesNumber){
    //setUpStatesCartesian(EnergyCutOff);
    //setUpStatesPolar(EnergyCutOff, h_omega, ParticlesNumber);
    setUpStatesPolarSorted(EnergyCutOff, h_omega, ParticlesNumber);

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

void QuantumDot::setUpStatesPolar(int EnergyCutOff, double h_omega, int ParticlesNumber) {
    //QuantumState * m_q_state = new QuantumState;
    homega = h_omega;
    NumberOfParticles = ParticlesNumber;
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

void QuantumDot::setUpStatesPolarSorted(int EnergyCutOff, double h_omega, int ParticlesNumber) {
    //QuantumState * m_q_state = new QuantumState;
    homega = h_omega;
    NumberOfParticles = ParticlesNumber;
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
           //m_q_state.flipSpin();
           //m_shells.push_back(m_q_state);
       }
       //negative m
       for (int j=1; j<=EnergyCutOff - oddShells[index]; j++){
           //cout << "[ " << index << " , " << j << " ]" << endl;
           n = index;
           m = -1*j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
           //m_q_state.flipSpin();
           //m_shells.push_back(m_q_state);
       }
    }

    //now sort m_shells vector
    pair<int,double> mapping;
    vector<pair<int,double>> vector_to_sort;
    for(int i = 0; i < m_shells.size(); i++) {
        QuantumState quantum_state = m_shells.at(i);
        double EnergyOfState = homega*((double)2.0*quantum_state.n() + abs(quantum_state.m()) + 1.0);
        mapping = make_pair(i, EnergyOfState);
        vector_to_sort.push_back(mapping);
        //cout << i << "   " << EnergyOfState << endl;
        //cout << vector_to_sort[i] << endl;
    }

    sort(vector_to_sort.begin(),vector_to_sort.end(),myComparison);
    vector<QuantumState> sorted_states;
    for(int i = 0; i < vector_to_sort.size(); i++) {
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
        m_shells.at(vector_to_sort.at(i).first).flipSpin();
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
        //cout << vector_to_sort.at(i).second << "    "  << vector_to_sort.at(i).first << endl;
    }
    m_shells = sorted_states;

}


void QuantumDot::setCoefficientMatrix(arma::mat CoefficientMatrix){
     m_C = CoefficientMatrix;
}

arma::mat QuantumDot::computeDensityMatrix(){
    int NumberOfStates = m_shells.size();
    //arma::mat C(NumberOfStates, NumberOfStates);
    //C.eye();
    arma::mat DensityMatrix(NumberOfStates, NumberOfStates);
    double sum = 0.0;
    for(int k = 0; k < NumberOfStates; k++) {
        for(int l = 0; l < NumberOfStates; l++) {
            for (int i=0; i < NumberOfParticles; i++) {
                    sum += m_C(k,i)*m_C(l,i);
                    DensityMatrix(k, l) = sum;
            }
        }
    }
    return DensityMatrix;
}



double QuantumDot::CalculateNonIntEnergy(int n, int m){
    return (2.0*(double)n + (double)abs(m) + 1.0)*homega;
}

void QuantumDot::computeHFmatrix(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    //arma::mat HF(NumberOfStates, NumberOfStates);
    m_HF.zeros(NumberOfStates,NumberOfStates);
    double FockElement;
    //arma::mat C(NumberOfStates, NumberOfStates);
    //C.eye();
    //setCoefficientMatrix(C);
    //arma::mat DensityMatrix = computeDensityMatrix();


    /*arma::mat C(NumberOfStates, NumberOfStates);
    C.eye();
    arma::mat DensityMatrix(NumberOfStates, NumberOfStates);
    DensityMatrix.zeros();

    double sum = 0.0;
    for(int k = 0; k < NumberOfStates; k++) {
    //QuantumState quantum_state_gama1 = m_shells.at(k);
    //int gama_n = quantum_state_gama1.n();
    //int gama_m = quantum_state_gama1.m();
    //int gama_sm = quantum_state_gama1.sm();
    //int gama_s = quantum_state_gama1.s();
        for(int l = 0; l < NumberOfStates; l++) {
            //QuantumState quantum_state_delta1 = m_shells.at(l);
            //int delta_n = quantum_state_delta1.n();
            //int delta_m = quantum_state_delta1.m();
            //int delta_sm = quantum_state_delta1.sm();
            //int delta_s = quantum_state_delta1.s();
                for (int i=0; i <NumberOfParticles; i++){
                    sum += C(k,i)*C(l,i);
                    DensityMatrix(k, l) = sum;
                }
            }
      }
*/
    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.sm();
        //int alpha_s = quantum_state_alpha.s();
        double nonInteractingPart = CalculateNonIntEnergy( alpha_n, alpha_m);
        m_HF(i, i) = nonInteractingPart;
        for(int j = 0; j < NumberOfStates; j++) {
            QuantumState quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            //int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.sm();
            //int beta_s = quantum_state_beta.s();
            for(int k = 0; k < NumberOfStates; k++) {
                QuantumState quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                //int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.sm();
                //int gama_s = quantum_state_gama.s();
                for(int l = 0; l < NumberOfStates; l++) {
                    QuantumState quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    //int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.sm();
                    //int delta_s = quantum_state_delta.s();
                    if (alpha_sm == beta_sm && gama_sm == delta_sm){
                        double TBME = Coulomb_HO(homega, alpha_n, alpha_sm, beta_n, beta_sm, gama_n, gama_sm, delta_n, delta_sm);
                        //cout << setprecision(12);
                        //cout << "< " << alpha_n << "," << alpha_sm << " ; " << beta_n << "," << gama_sm << " || V || " << gama_n << "," << gama_sm << " ; " << delta_n << "," << delta_sm << " > = " << TBME << std::endl;
                        //cout << "< " << p_n << "," << p_sm << " ; " << q_n << "," << q_sm << " || V || " << r_n << "," << r_sm << " ; " << s_n << "," << s_sm << " > = " << endl;
                        //HF(index_1,index_2) += density()*TBME;

                        FockElement += DensityMatrix(k,l)*TBME;
                        m_HF(i, j) += FockElement;
                    } else {
                        m_HF(i, j) = 0.0;
                    }
                }
            }
        }
    }
    //m_HF.print();
}


double QuantumDot::computeHartreeFockEnergyDifference(){
    return ((arma::accu(abs(eigval - eigval_previous)))/(double)m_shells.size());
    //return energy_diff;
}

void QuantumDot::computeHartreeFockEnergy(){
    //cout << (arma::accu((eigval))/(double)m_shells.size()) << endl;
    int NumberOfStates = m_shells.size();
    int FermiLevel = 3;
    //m_HF.zeros(NumberOfStates,NumberOfStates);
    double selfConsistentFIeldIterations;
    double ExchangePart;

    double SingleParticleEnergies = 0.0;
    for(int f = 0; f < FermiLevel; f++){
        SingleParticleEnergies += eigval(f);
    }


    for(int a = 0; a < FermiLevel; a++) {
        for(int b = 0; b < FermiLevel; b++) {
            for(int i = 0; i < NumberOfParticles; i++) {
                QuantumState quantum_state_alpha = m_shells.at(i);
                int alpha_n = quantum_state_alpha.n();
                //int alpha_m = quantum_state_alpha.m();
                int alpha_sm = quantum_state_alpha.sm();
                for(int j = 0; j < NumberOfParticles; j++) {
                    QuantumState quantum_state_beta = m_shells.at(j);
                    int beta_n = quantum_state_beta.n();
                    //int beta_m = quantum_state_beta.m();
                    int beta_sm = quantum_state_beta.sm();
                    //int beta_s = quantum_state_beta.s();
                    for(int k = 0; k < NumberOfParticles; k++) {
                        QuantumState quantum_state_gama = m_shells.at(k);
                        int gama_n = quantum_state_gama.n();
                        //int gama_m = quantum_state_gama.m();
                        int gama_sm = quantum_state_gama.sm();
                        //int gama_s = quantum_state_gama.s();
                        for(int l = 0; l < NumberOfParticles; l++) {
                            QuantumState quantum_state_delta = m_shells.at(l);
                            int delta_n = quantum_state_delta.n();
                            //int delta_m = quantum_state_delta.m();
                            int delta_sm = quantum_state_delta.sm();
                            //int delta_s = quantum_state_delta.s();
                            if (alpha_sm == beta_sm && gama_sm == delta_sm){
                                double TBME = Coulomb_HO(homega, alpha_n, alpha_sm, beta_n, beta_sm, gama_n, gama_sm, delta_n, delta_sm);
                                selfConsistentFIeldIterations += m_C(a,i)*m_C(a,j)*m_C(b,k)*m_C(b,l)*TBME;
                                ExchangePart += selfConsistentFIeldIterations;
                            } else {
                                ExchangePart += 0.0;
                            }
                        }
                    }
                }
            }
        }
    }



    cout << "SPEnergies" << SingleParticleEnergies << endl;
    cout << "Exchange" << ExchangePart << endl;
    cout << "Total" << SingleParticleEnergies - 0.5*ExchangePart << endl;
}


void QuantumDot::applyHartreeFockMethod(){
    int NumberOfStates = m_shells.size();
    arma::mat C(NumberOfStates, NumberOfStates);

    C.eye();
    setCoefficientMatrix(C);
    double difference = 10; //dummy value to handle first iteration
    double epsilon = 10e-8;

    //arma::cx_vec eigval;
    //arma::cx_mat eigvec;
    eigval_previous.zeros(NumberOfStates);
    int i = 0;
    while (epsilon < difference && i < 10000){
        arma::mat x_DensityMatrix = computeDensityMatrix();
        computeHFmatrix(x_DensityMatrix);
        arma::eig_sym(eigval, eigvec, m_HF);
        //eigval.print();
        setCoefficientMatrix(eigvec);
        difference = computeHartreeFockEnergyDifference();
        eigval_previous = eigval;
        i++;

    }
    eigval.print();
    eigvec.print();
    computeHartreeFockEnergy();
}


void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << "n = " << quantum_state.n() << endl;
        cout << "m = " <<quantum_state.m() << endl;
        cout << "sm = " <<quantum_state.sm() << endl;
        cout << "s = " <<quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

void QuantumDot::getQuantumDotStatesNumber(){
    cout << "Number of available states of system is " << m_shells.size() << endl;
}

