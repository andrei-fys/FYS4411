#include "quantumdot.h"
#include "quantumforce.h"
#include <iostream>
#include <cmath>
#include <iomanip>  //on Mac setprecision
#include <fstream>

using namespace std;



QuantumDot::QuantumDot(double h_omega, int ParticlesNumber){
    initialize(ParticlesNumber);
    homega = h_omega;
}

void QuantumDot::setVariationalParameters(double first, double second){
    alpha = first;
    beta = second;
}


void QuantumDot::initialize(int ParticlesNumber){
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    for(int i=0; i<ParticlesNumber; i++) {
        Particle *particle = new Particle();
        double x = RandomNumberGenerator(gen); // random number in the interval [0,1]
        double y = RandomNumberGenerator(gen);
        double z = 0;
        particle->position.set(x,y,z);
        m_particles.push_back(particle);
    }
}

void QuantumDot::applyVMC(int MCSamples){
    string ResultsFile = "VMC_ImpSampl";
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    normal_distribution<double> GaussianBlur(0.0,1);
    m_localEnergy.clear();
    int MC_counter = 0;
    double MeanLocalEnergySquared = 0;
    double MeanLocalEnergy = 0;
    int accept=0;
    for(int i=0; i<MCSamples; i++){
        for(int j=0; j< m_particles.size() ; j++) {
            Particle *particle = m_particles[j];
            Particle *particlemove;
            for(int k=0; k< m_particles.size() ; k++) {
                if (k != j){
                    particlemove = m_particles[k];
                }
            }

            QuantumForce QForce(alpha, homega);
            double QForceOld = QForce.evaluate(particle->position[0],particle->position[1],
                    particlemove->position[0],particlemove->position[1]);
            double addition1 =  GaussianBlur(gen)*sqrt(dt) +QForceOld*dt*D;
            double addition2 =  GaussianBlur(gen)*sqrt(dt) +QForceOld*dt*D;

            double X_new = particle->position.x() + addition1;
            double Y_new = particle->position.y() + addition2;


            particle->positionNew.setX(X_new);
            particle->positionNew.setY(Y_new);
            //start of Metropolis-Hasting´s
            double Rold2 = particle->position.lengthSquared();
            double Rnew2 = particle->positionNew.lengthSquared();
            double w = exp(alpha*homega*(Rold2 - Rnew2));
            w *= calculateTransitionProbabilityImpSampl(j);
            double r = RandomNumberGenerator(gen);
            if (w >= r) {
               particle->position.setX(X_new);
               particle->position.setY(Y_new);
               accept++;
            }
        }
        MC_counter++;
        double LocalEnergy = calculateLocalEnergy();
        //writeLocalEnergyToFile(LocalEnergy, ResultsFile);
        m_localEnergy.push_back(LocalEnergy);
        MeanLocalEnergySquared += LocalEnergy*LocalEnergy;
        MeanLocalEnergy += LocalEnergy;
    }
    double Energy =  (double)MeanLocalEnergy/MCSamples;
    double Energy2 = (double)MeanLocalEnergySquared/MCSamples;
    cout << "Energy "  << Energy << endl;
    cout << "Energy2 " << Energy2 << endl;
    cout << "Variance " << (Energy2 - Energy*Energy)/MCSamples << endl;
    cout << "Accept " << accept/2 << endl;
    cout << "MCcounter " << MC_counter << endl;
    writeVectorToFile(ResultsFile);
}

void QuantumDot::applyVMCstandard(int MCSamples){
    string ResultsFile = "VMC_Standatd";
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
    uniform_real_distribution<double> RandomNumberGenerator1(0.0,1.0);
    m_localEnergy.clear();
    int MC_counter = 0;
    double h = 1.3;
    //double h = 0.005;
    double MeanLocalEnergySquared = 0;
    double MeanLocalEnergy = 0;
    int accept=0;
    for(int i=0; i<MCSamples; i++){ // MC sampling
        for(Particle *particle : m_particles){
            double X_new = particle->position.x() + (RandomNumberGenerator(gen)*h);
            double Y_new = particle->position.y() + (RandomNumberGenerator(gen)*h);
            particle->positionNew.setX(X_new);
            particle->positionNew.setY(Y_new);
            //start of Metropolis
            double Rold2 = particle->position.lengthSquared();
            double Rnew2 = particle->positionNew.lengthSquared();
            double w = exp(alpha*homega*(Rold2 - Rnew2));
            double r = RandomNumberGenerator1(gen);
            if (w >= r) {
               particle->position.setX(X_new);
               particle->position.setY(Y_new);
               accept++;
            }
        }
        MC_counter++;
        double LocalEnergy = calculateLocalEnergy();
        m_localEnergy.push_back(LocalEnergy);
        //writeLocalEnergyToFile(LocalEnergy, ResultsFile);
        MeanLocalEnergySquared += LocalEnergy*LocalEnergy;
        MeanLocalEnergy += LocalEnergy;
    }  //End of MC sampling
    double Energy =  (double)MeanLocalEnergy/MCSamples;
    double Energy2 = (double)MeanLocalEnergySquared/MCSamples;
    cout << "Energy "  << Energy << endl;
    cout << "Energy2 " << Energy2 << endl;
    cout << "Variance " << (Energy2 - Energy*Energy)/MCSamples << endl;
    cout << "Accept " << accept/2 << endl;
    cout << "MCcounter " << MC_counter << endl;
    writeVectorToFile(ResultsFile);
}

double QuantumDot::calculateTransitionProbabilityImpSampl(int j){
    double QForceOld;
    double QForceNew;
    double LengthOld = 0;
    double LengthNew = 0;
    QuantumForce QForce(alpha, homega);
    for(int i=0; i< m_particles.size() ; i++) {
        Particle *particle = m_particles[i];
        Particle *moving_particle = m_particles[j];
        if (i != j){
            QForceNew = QForce.evaluate(particle->position[0],particle->position[1],
                    moving_particle->positionNew[0],moving_particle->positionNew[1]);
            QForceOld = QForce.evaluate(particle->position[0],particle->position[1],
                          moving_particle->position[0],moving_particle->position[1]);
            LengthOld = moving_particle->position.length();
            LengthNew = moving_particle->positionNew.length();
        }
    }


    double GreenPart = (0.25*D*dt*(QForceOld*QForceOld - QForceNew*QForceNew)) + 0.5*(LengthOld - LengthNew)*(QForceOld + QForceNew);
    return exp(GreenPart);
}





/*double computeRelativeDistance(){
    for(int i=0; i<QuantumDot.m_particles.size(); i++) {
       Particle *part_i = QuantumDot.m_particles[i];
       for(int j=i+1; j<QuantumDot.m_particles.size(); j++) {
           Particle *part_j = QuantumDot.m_particles[j];

           double dx = part_i->position[0] - part_j->position[0];
           double dy = part_i->position[1] - part_j->position[1];
           return sqrt(dx*dx + dy*dy);
       }
    }
}
*/

double QuantumDot::calculateLocalEnergy(){
    double R2 = 0;
    vec3 RelativeDistance;
    if (m_Coulomb == 0) {
        for(Particle *particle : m_particles){
            double length2 = particle->position.lengthSquared();
            R2 += length2;
        }
        return 0.5*homega*homega*(1.0 - alpha*alpha)*R2 +2.0*alpha*homega;
    } else {
        for(int i=0; i<m_particles.size(); i++) {
            Particle *particle_i = m_particles[i];
            double length2 = particle_i->position.lengthSquared();
            R2 += length2;
            for(int j=0; j<m_particles.size(); j++) {
                Particle *particle_j = m_particles[j];
                if (i < j){
                    RelativeDistance = particle_i->position - particle_j->position;
                }
            }
        }
    return 0.5*homega*homega*(1.0 - alpha*alpha)*R2 +2.0*alpha*homega + 1.0/RelativeDistance.length();
    }
}



void QuantumDot::getQuantumDotParticlesCoordinates(){
    for(Particle *particle : m_particles){
        particle->position.print();
    }
}

void QuantumDot::setUpStatesCartesian(int ParticlesNumber) {
    QuantumState m_q_state;
    int nx, ny;
    for (int i = 1; i < ParticlesNumber + 1; i++){
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

void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << "n = " << quantum_state.nx() << endl;
        cout << "m = " <<quantum_state.ny() << endl;
        cout << "sm = " <<quantum_state.sm() << endl;
        cout << "s = " <<quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

void QuantumDot::getQuantumDotStatesNumber(){
    cout << "Number of electrons is " << m_shells.size() << endl;
}


void QuantumDot::setCoulombInterraction(int trigger){
    m_Coulomb = trigger;
}

void QuantumDot::writeLocalEnergyToFile(double LocalEnergy, string ResultsFile){
    ofstream ofile;
    ofile.open(ResultsFile, ios::app);
    ofile << setprecision(12);
    ofile << LocalEnergy << endl;
    ofile.close();
}

void QuantumDot::writeVectorToFile(string ResultsFile){
    ofstream ofile;
    ofile.open(ResultsFile, ios::app);
    ofile << setprecision(12);
    for(double element : m_localEnergy){
        ofile << element << endl;
    }
    ofile.close();
}
