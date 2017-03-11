#include "quantumdot.h"
#include <iostream>


using namespace std;

bool myComparison(const pair<int,double> &a,const pair<int,double> &b)
{
       return a.second<b.second;
}

QuantumDot::QuantumDot(int EnergyCutOff, double h_omega, int ParticlesNumber){
    setUpStatesPolarSorted(EnergyCutOff, h_omega, ParticlesNumber);
}


void QuantumDot::setUpStatesCartesian(int EnergyCutOff) {
    QuantumState m_q_state;
    int nx, ny;
    for (int i = 1; i < EnergyCutOff + 1; i++){
        for (int j = 0; j < i; j++){
            nx = j;
            ny = i - nx -1;
            m_q_state.set(nx, ny, m_sm, m_s);
            m_shells.push_back(m_q_state);
            m_q_state.flipSpin();
            m_shells.push_back(m_q_state);
         }
    }
}

void QuantumDot::setUpStatesPolarSorted(int EnergyCutOff, double h_omega, int ParticlesNumber) {
    homega = h_omega;
    NumberOfParticles = ParticlesNumber;
    m_EnergyCutOff = EnergyCutOff;
    QuantumState m_q_state;
    vector<int> oddShells;
    int n, m;

    //loop to find all odd shell numbers including energyCutOff shell
    for (int i = 1; i <= EnergyCutOff; i++){
        if (i % 2 != 0){
            oddShells.push_back(i);
        }
    }
    for (int index = 0; index < oddShells.size(); ++index){
       //positive m
       for (int j=0; j<=EnergyCutOff - oddShells[index]; j++){
           n = index;
           m = j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
       }
       //negative m
       for (int j=1; j<=EnergyCutOff - oddShells[index]; j++){
           n = index;
           m = -1*j;
           m_q_state.set(n, m, m_sm, m_s);
           m_shells.push_back(m_q_state);
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
    }

    sort(vector_to_sort.begin(),vector_to_sort.end(),myComparison);
    vector<QuantumState> sorted_states;
    for(int i = 0; i < vector_to_sort.size(); i++) {
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
        m_shells.at(vector_to_sort.at(i).first).flipSpin();
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
    }
    m_shells = sorted_states;
    CalculateNonIntEnergy();

}


void QuantumDot::setCoefficientMatrix(arma::mat CoefficientMatrix){
     m_C = CoefficientMatrix;
}

arma::mat QuantumDot::computeDensityMatrix(){
    int NumberOfStates = m_shells.size();
    arma::mat DensityMatrix(NumberOfStates, NumberOfStates);

    for(int k = 0; k < NumberOfStates; k++) {
        for(int l = 0; l < NumberOfStates; l++) {
            double sum = 0.0;
            for (int i=0; i < NumberOfParticles; i++) {
                    sum += m_C(k,i)*m_C(l,i);
                    DensityMatrix(k, l) = sum;
            }
        }
    }
    return DensityMatrix;
}

void QuantumDot::CalculateNonIntEnergy(){
    int NumberOfStates = m_shells.size();
    m_HOEnergies.zeros(NumberOfStates,NumberOfStates);
    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state = m_shells.at(i);
        m_HOEnergies(i, i) = (2.0*(double)quantum_state.n() + (double)abs(quantum_state.m()) + 1.0)*homega;
    }
}

void QuantumDot::computeHFmatrix(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    m_HF.zeros(NumberOfStates,NumberOfStates);
    double FockElement = 0;

    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.sm();

        for(int j = 0; j < NumberOfStates; j++) {
            QuantumState quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.sm();

            for(int k = 0; k < NumberOfStates; k++) {
                QuantumState quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.sm();

                for(int l = 0; l < NumberOfStates; l++) {
                    QuantumState quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.sm();
                    double TBME = 0.0;
                    double tbme1 = 0.0;
                    double tbme2 = 0.0;

                    if ((alpha_sm == beta_sm && gama_sm == delta_sm)){ /*&&
                            (alpha_sm + gama_sm == beta_sm + delta_sm) &&
                                (alpha_m + gama_m == beta_m + delta_m)){*/

                        tbme1 = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, beta_n, beta_m,  delta_n, delta_m);
                    }
                    if ((alpha_sm == delta_sm && gama_sm == beta_sm)){ /*&&
                            (alpha_sm + gama_sm == beta_sm + delta_sm) &&
                                (alpha_m + gama_m == delta_m + beta_m)){*/
                        tbme2   =  Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, delta_n, delta_m, beta_n, beta_m);
                    }
                    TBME = tbme1 - tbme2;
                    FockElement += DensityMatrix(k,l)*TBME;
                }
            }
            if (i == j) {
                m_HF(i, i) += m_HOEnergies(i, i);
            }
            m_HF(i, j) += FockElement;
            FockElement = 0.0;
        }
    }
}


double QuantumDot::computeHartreeFockEnergyDifference(){
    return ((arma::accu(abs(eigval - eigval_previous)))/(double)m_shells.size());
}

void QuantumDot::computeHartreeFockEnergy(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    int FermiLevel = NumberOfParticles;
    double selfConsistentFIeldIterations = 0.0;
    double ExchangePart = 0.0;
    double SingleParticleEnergies = 0.0;

    for(int f = 0; f < FermiLevel; f++){
        SingleParticleEnergies += eigval(f);
    }

    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.sm();

        for(int j = 0; j < NumberOfStates; j++) {
            QuantumState quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.sm();

            for(int k = 0; k < NumberOfStates; k++) {
                QuantumState quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.sm();

                for(int l = 0; l < NumberOfStates; l++) {
                    QuantumState quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.sm();

                    double TBME = 0.0;
                    double tbme1 = 0.0;
                    double tbme2 = 0.0;
                    if ((alpha_sm == beta_sm) && (gama_sm == delta_sm)){
                       tbme1 = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, beta_n, beta_m,  delta_n, delta_m);
                    }
                    if ((alpha_sm == delta_sm) && (gama_sm == beta_sm)){
                       tbme2   =  Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, delta_n, delta_m, beta_n, beta_m);
                    }
                    TBME = tbme1 - tbme2;
                    selfConsistentFIeldIterations = DensityMatrix(i,j)*DensityMatrix(k,l)*TBME;
                    ExchangePart += selfConsistentFIeldIterations;
                }
            }
        }
    }
    double HF_Energy = SingleParticleEnergies - 0.5*ExchangePart;
    // Uncoment for debug
    //cout << "SPEnergies " << SingleParticleEnergies << endl;
    //cout << "Exchange " << ExchangePart << endl;
    cout << "===================================================================" << endl;
    cout << "Num of electrons = " << NumberOfParticles << endl;
    cout << "Num of shells = " << m_EnergyCutOff << endl;
    cout << "Omega = " << homega << endl;
    cout << "Total energy " << HF_Energy << endl;
    writeToFile(HF_Energy, NumberOfParticles, m_EnergyCutOff, homega);
}

void QuantumDot::applyHartreeFockMethod(){
    int NumberOfStates = m_shells.size();
    arma::mat C(NumberOfStates, NumberOfStates);

    C.eye();
    setCoefficientMatrix(C);
    double difference = 10; //dummy value to handle first iteration
    double epsilon = 10e-10;

    eigval_previous.zeros(NumberOfStates);
    int i = 0;
    while (epsilon < difference && i < 1000){
        arma::mat x_DensityMatrix = computeDensityMatrix();
        computeHFmatrix(x_DensityMatrix);
        arma::eig_sym(eigval, eigvec, m_HF);
        setCoefficientMatrix(eigvec);
        difference = computeHartreeFockEnergyDifference();
        eigval_previous = eigval;
        i++;

    }

    arma::mat y_DensityMatrix = computeDensityMatrix();
    computeHartreeFockEnergy(y_DensityMatrix);
    cout << "Number of iterations " << i << endl;
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

void QuantumDot::writeToFile(double HF_Energy, int NumberOfParticles, int m_EnergyCutOff, double homega){
    ofstream ofile;
    ofile.open(ResultsFile, ios::app);
    ofile << "===============================" << endl;
    ofile << "Num of electrons = " << NumberOfParticles << endl;
    ofile << "Num of shells = " << m_EnergyCutOff << endl;
    ofile << "Omega = " << homega << endl;
    ofile << "Total energy " << HF_Energy << endl;
    ofile << "Eigenvalues: " << endl;
    for (int i = 0; i < eigval.size(); i++){
        ofile << eigval(i) << "     "<< m_HOEnergies(i, i) << endl;
    }
    ofile.close();
}
