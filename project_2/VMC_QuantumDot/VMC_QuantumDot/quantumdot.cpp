#include "quantumdot.h"
//#include "quantumforce.h"
#include <iostream>
#include <cmath>
#include <iomanip>  //on Mac setprecision
#include <fstream>
//#include "mpi.h"

using namespace std;

QuantumDot::QuantumDot(double h_omega, int ParticlesNumber){
    initialize(ParticlesNumber);
    homega = h_omega;
    m_homega2 = h_omega*h_omega;
    m_InsideSteepestDescent = false;
    m_ExpectationLocalEnergyDerivativeAlphaSecondTerm = 0.0;
    m_ExpectationLocalEnergyDerivativeBetaSecondTerm = 0.0;
    m_ExpectationLocalEnergyDerivativeAlphaFirstTerm = 0.0;
    m_ExpectationLocalEnergyDerivativeBetaFirstTerm = 0.0;
    m_LocalEnergyWFDerivativeAlpha = 0.0;
    m_LocalEnergyWFDerivativeBeta = 0.0;
    m_MeanLocalEnergyWFDerivativeAlpha = 0.0;
    m_MeanLocalEnergyWFDerivativeBeta = 0.0;
    m_ExpectationLocalEnergyDerivativeAlpha = 0.0;
    m_ExpectationLocalEnergyDerivativeBeta = 0.0;
    //m_Energy = 0.0;
}

void QuantumDot::setVariationalParameters(double first, double second){
    alpha = first;
    m_alpha2 = first*first;
    beta = second;
    m_alphaomega = alpha*homega;
    m_alphaomega2 = m_alphaomega*m_alphaomega;
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

/*void QuantumDot::applyVMCMPI(int MCSamples){
    int numprocs, my_rank;
    numprocs = 8;
    //MPI Part
    MPI_Init (NULL, NULL);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    //MPI Part
    string ResultsFile = "VMC_ImpSampl";
    random_device rd;
    mt19937_64 gen(rd());
    gen.discard(700000);
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    normal_distribution<double> GaussianBlur(0.0, 1.0);
    m_localEnergy.clear();
    int MC_counter = 0;
    double MeanLocalEnergySquared = 0;
    double MeanLocalEnergy = 0;
    int accept=0;
    //cout << "alpha "  << alpha << endl;
    //cout << "beta " << beta << endl;
    for(int i=0; i<MCSamples; i++){
        for(size_t j=0; j< m_particles.size() ; j++) {
            Particle *particle = m_particles[j];
            Particle *particlemove;
            for(size_t k=0; k< m_particles.size() ; k++) {
                if (k != j){
                    particlemove = m_particles[k];
                }
            }

            //QuantumForce QForce(alpha, homega);
            //double QForceOld = QForce.evaluate(particle->position[0],particle->position[1],
            //        particlemove->position[0],particlemove->position[1]);
            calculateQuantumForce(j);
            double additionX =  GaussianBlur(gen)*sqrt(dt) +m_QForceOld[0]*dt*D;
            double additionY =  GaussianBlur(gen)*sqrt(dt) +m_QForceOld[1]*dt*D;

            double X_new = particle->position.x() + additionX;
            double Y_new = particle->position.y() + additionY;


            particle->positionNew.setX(X_new);
            particle->positionNew.setY(Y_new);
            calculateQuantumForceNew(j);
            //start of Metropolis-Hasting´s
            double Rold2 = particle->position.lengthSquared();
            double Rnew2 = particle->positionNew.lengthSquared();
            double w = exp(alpha*homega*(Rold2 - Rnew2))*calculateJastrowRatio(j);
            w *= calculateGreenFunctionRatio(j);
            double r = RandomNumberGenerator(gen);
            if (w >= r) {
               particle->position.setX(X_new);
               particle->position.setY(Y_new);
               accept++;
            }
        }
        MC_counter++;
        double LocalEnergy = calculateLocalEnergy();
        if (m_InsideSteepestDescent = true){
            m_MeanLocalEnergyWFDerivativeAlpha += m_LocalEnergyWFDerivativeAlpha;
            m_MeanLocalEnergyWFDerivativeBeta += m_LocalEnergyWFDerivativeBeta;
            m_ExpectationLocalEnergyDerivativeAlphaFirstTerm += m_LocalEnergyWFDerivativeAlpha*LocalEnergy;
            m_ExpectationLocalEnergyDerivativeBetaFirstTerm += m_LocalEnergyWFDerivativeBeta*LocalEnergy;
        }
        //writeLocalEnergyToFile(LocalEnergy, ResultsFile);
        //m_localEnergy.push_back(LocalEnergy);
        MeanLocalEnergySquared += LocalEnergy*LocalEnergy;
        MeanLocalEnergy += LocalEnergy;

    }
    double Energy =  (double)MeanLocalEnergy/MCSamples;
    double Energy2 = (double)MeanLocalEnergySquared/MCSamples;
    double accept_slave = ((double)accept/(double)(2.0*MCSamples))*100.0;
    double variance_slave = (Energy2 - Energy*Energy)/MCSamples;
    if (m_InsideSteepestDescent = true){
        m_ExpectationLocalEnergyDerivativeAlphaSecondTerm = Energy*((double)m_MeanLocalEnergyWFDerivativeAlpha/MCSamples);
        m_ExpectationLocalEnergyDerivativeBetaSecondTerm = Energy*((double)m_MeanLocalEnergyWFDerivativeBeta/MCSamples);
        m_ExpectationLocalEnergyDerivativeAlphaFirstTerm /= (double)MCSamples;
        m_ExpectationLocalEnergyDerivativeBetaFirstTerm /= (double)MCSamples;
        m_ExpectationLocalEnergyDerivativeAlpha = 2.0*(m_ExpectationLocalEnergyDerivativeAlphaFirstTerm-m_ExpectationLocalEnergyDerivativeAlphaSecondTerm);
        m_ExpectationLocalEnergyDerivativeBeta = 2.0*(m_ExpectationLocalEnergyDerivativeBetaFirstTerm-m_ExpectationLocalEnergyDerivativeBetaSecondTerm);
    }
    //m_Energy = Energy;
    // ** cout << "Energy "  << Energy << endl;
    //cout << "Energy2 " << Energy2 << endl;
    //cout << "Variance " << (Energy2 - Energy*Energy)/MCSamples << endl;
    // ** cout << "Accept " << ((double)accept/(double)(2.0*MCSamples))*100.0 << endl;
    //cout << "MCcounter " << MC_counter << endl;
    //writeVectorToFile(ResultsFile);
    double master_Energy = 0.0;
    double accept_master = 0.0;
    double variance_master = 0.0;
    MPI_Reduce(&Energy, &master_Energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accept_slave, &accept_master, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&variance_slave, &variance_master, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_rank == 0){
        cout << "alpha "  << alpha << endl;
        cout << "beta " << beta << endl;
        cout << "Energy " << (double) master_Energy/numprocs << endl;
        cout << "Variance " << (double) variance_master/numprocs << endl;
        cout << "Accept % :" << (double) accept_master/numprocs << endl;
        cout << "Total number of MC samples " << MCSamples*numprocs << endl;

    }
    MPI_Finalize();
}
*/
void QuantumDot::applyVMC(int MCSamples){
    string ResultsFile = "VMC_ImpSampl";
    random_device rd;
    mt19937_64 gen(rd());
    gen.discard(700000);
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    normal_distribution<double> GaussianBlur(0.0, 1.0);
    m_localEnergy.clear();
    int MC_counter = 0;
    double MeanLocalEnergySquared = 0;
    double MeanLocalEnergy = 0;
    int accept=0;
    //cout << "alpha "  << alpha << endl;
    //cout << "beta " << beta << endl;
    for(int i=0; i<MCSamples; i++){
        for(size_t j=0; j< m_particles.size() ; j++) {
            Particle *particle = m_particles[j];
            Particle *particlemove;
            for(size_t k=0; k< m_particles.size() ; k++) {
                if (k != j){
                    particlemove = m_particles[k];
                }
            }

            //QuantumForce QForce(alpha, homega);
            //double QForceOld = QForce.evaluate(particle->position[0],particle->position[1],
            //        particlemove->position[0],particlemove->position[1]);
            calculateQuantumForce(j);
            double additionX =  GaussianBlur(gen)*sqrt(dt) +m_QForceOld[0]*dt*D;
            double additionY =  GaussianBlur(gen)*sqrt(dt) +m_QForceOld[1]*dt*D;

            double X_new = particle->position.x() + additionX;
            double Y_new = particle->position.y() + additionY;


            particle->positionNew.setX(X_new);
            particle->positionNew.setY(Y_new);
            calculateQuantumForceNew(j);
            //start of Metropolis-Hasting´s
            double Rold2 = particle->position.lengthSquared();
            double Rnew2 = particle->positionNew.lengthSquared();
            double w = exp(alpha*homega*(Rold2 - Rnew2))*calculateJastrowRatio(j);
            w *= calculateGreenFunctionRatio(j);
            double r = RandomNumberGenerator(gen);
            if (w >= r) {
               particle->position.setX(X_new);
               particle->position.setY(Y_new);
               accept++;
            }
        }
        MC_counter++;
        double LocalEnergy = calculateLocalEnergy();
        if (m_InsideSteepestDescent = true){
            m_MeanLocalEnergyWFDerivativeAlpha += m_LocalEnergyWFDerivativeAlpha;
            m_MeanLocalEnergyWFDerivativeBeta += m_LocalEnergyWFDerivativeBeta;
            m_ExpectationLocalEnergyDerivativeAlphaFirstTerm += m_LocalEnergyWFDerivativeAlpha*LocalEnergy;
            m_ExpectationLocalEnergyDerivativeBetaFirstTerm += m_LocalEnergyWFDerivativeBeta*LocalEnergy;
        }
        //writeLocalEnergyToFile(LocalEnergy, ResultsFile);
        //m_localEnergy.push_back(LocalEnergy);
        MeanLocalEnergySquared += LocalEnergy*LocalEnergy;
        MeanLocalEnergy += LocalEnergy;

    }
    double Energy =  (double)MeanLocalEnergy/MCSamples;
    double Energy2 = (double)MeanLocalEnergySquared/MCSamples;
    if (m_InsideSteepestDescent = true){
        m_ExpectationLocalEnergyDerivativeAlphaSecondTerm = Energy*((double)m_MeanLocalEnergyWFDerivativeAlpha/MCSamples);
        m_ExpectationLocalEnergyDerivativeBetaSecondTerm = Energy*((double)m_MeanLocalEnergyWFDerivativeBeta/MCSamples);
        m_ExpectationLocalEnergyDerivativeAlphaFirstTerm /= (double)MCSamples;
        m_ExpectationLocalEnergyDerivativeBetaFirstTerm /= (double)MCSamples;
        m_ExpectationLocalEnergyDerivativeAlpha = 2.0*(m_ExpectationLocalEnergyDerivativeAlphaFirstTerm-m_ExpectationLocalEnergyDerivativeAlphaSecondTerm);
        m_ExpectationLocalEnergyDerivativeBeta = 2.0*(m_ExpectationLocalEnergyDerivativeBetaFirstTerm-m_ExpectationLocalEnergyDerivativeBetaSecondTerm);
    }
    //m_Energy = Energy;
    cout << "Alpha "  << alpha << endl;
    cout << "Beta "  << beta << endl;
    cout << "Energy "  << Energy << endl;
    //cout << "Energy2 " << Energy2 << endl;
    cout << "Variance " << (Energy2 - Energy*Energy)/MCSamples << endl;
    cout << "Accept " << ((double)accept/(double)(2.0*MCSamples))*100.0 << endl;
    cout << "MCcounter " << MC_counter << endl;
    //writeVectorToFile(ResultsFile);
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
    double MeanLocalEnergySquared = 0;
    double MeanLocalEnergy = 0;
    int accept=0;
    for(int i=0; i<MCSamples; i++){ // MC sampling
        //for(Particle *particle : m_particles){
        for(size_t j=0; j< m_particles.size() ; j++) {
            Particle *particle = m_particles[j];
            double X_new = particle->position.x() + (RandomNumberGenerator(gen)*h);
            double Y_new = particle->position.y() + (RandomNumberGenerator(gen)*h);
            particle->positionNew.setX(X_new);
            particle->positionNew.setY(Y_new);
            //start of Metropolis
            double Rold2 = particle->position.lengthSquared();
            double Rnew2 = particle->positionNew.lengthSquared();
            double w = exp(alpha*homega*(Rold2 - Rnew2));
            if (m_Jastrow == 1) {
                w *= calculateJastrowRatio(j);
            }
            double r = RandomNumberGenerator1(gen);
            if (w >= r) {
               particle->position.setX(X_new);
               particle->position.setY(Y_new);
               accept++;
            }
        }
        MC_counter++;
        double LocalEnergy = calculateLocalEnergy();
        //m_localEnergy.push_back(LocalEnergy);
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
    //writeVectorToFile(ResultsFile);
}

void QuantumDot::calculateQuantumForce(size_t j){
    vec3 RelativeDistance;
    double Fx = 0.0;
    double Fy = 0.0;
    double x2, y2;
    Particle *moving_particle = m_particles[j];
    for(size_t i=0; i< m_particles.size() ; i++) {
        Particle *particle = m_particles[i];
        if (i != j){
            RelativeDistance = particle->position - moving_particle->position;
            x2 = particle->position[0];
            y2 = particle->position[1];
        }
    }
    double R12 = RelativeDistance.length();
    double R12beta = (R12 + beta);
    double R12beta2 = 2.0*m_a/(R12beta*R12beta*R12);

    Fx = R12beta2*(moving_particle->position[0] - x2) - 2.0*m_alphaomega*moving_particle->position[0];
    Fy = R12beta2*(moving_particle->position[1] - y2) - 2.0*m_alphaomega*moving_particle->position[1];

    m_QForceOld.setX(Fx);
    m_QForceOld.setY(Fy);

}

void QuantumDot::calculateQuantumForceNew(size_t j){
    vec3 RelativeDistance;
    double x2, y2;
    Particle *moving_particle = m_particles[j];
    for(size_t i=0; i< m_particles.size() ; i++) {
        Particle *particle = m_particles[i];
        if (i != j){
            RelativeDistance = particle->position - moving_particle->positionNew;
            x2 = particle->position[0];
            y2 = particle->position[1];
        }
    }
    double R12 = RelativeDistance.length();
    double R12beta = (R12 + beta);
    double R12beta2 = 2.0*m_a/(R12beta*R12beta*R12);

    double Fx = R12beta2*(moving_particle->positionNew[0] - x2) - 2.0*m_alphaomega*moving_particle->positionNew[0];
    double Fy = R12beta2*(moving_particle->positionNew[1] - y2) - 2.0*m_alphaomega*moving_particle->positionNew[1];

    m_QForceNew.setX(Fx);
    m_QForceNew.setY(Fy);

}


double QuantumDot::calculateGreenFunctionRatio(size_t j){
    Particle *particle = m_particles[j];
    vec3 QForceDiff;
    QForceDiff = m_QForceOld - m_QForceNew;
    double GreenFun = D*dt*0.25*QForceDiff.lengthSquared() + (particle->position[0] - particle->positionNew[0])*(m_QForceOld[0] - m_QForceNew[0])*0.5 +
    (particle->position[1] - particle->positionNew[1])*(m_QForceOld[1] - m_QForceNew[1])*0.5;
    return exp(GreenFun);
}

double QuantumDot::calculateJastrowRatio(size_t j){
    if (m_Jastrow == 0) {
        return 1.0;
    } else {
        Particle *moving_particle = m_particles[j];
        vec3 RelativeDistance_old;
        vec3 RelativeDistance_new;
        for(size_t i=0; i< m_particles.size() ; i++) {
            Particle *particle = m_particles[i];
            if (i != j){
                RelativeDistance_old = particle->position - moving_particle->position;
                RelativeDistance_new = particle->position - moving_particle->positionNew;
            }
        }
        return exp(2.0*m_a*(1.0/(beta + 1.0/RelativeDistance_new.length()) - 1.0/(beta + 1.0/RelativeDistance_old.length())));
    }
}

double QuantumDot::calculateLocalEnergy(){
    double R2 = 0;
    vec3 RelativeDistance;
    if (m_Jastrow == 0) {
        calculateLocalEnergyWithoutJastrow();
    } else {
        for(size_t i=0; i<m_particles.size(); i++) {
            Particle *particle_i = m_particles[i];
            double length2 = particle_i->position.lengthSquared();
            R2 += length2;
            for(size_t j=0; j<m_particles.size(); j++) {
                Particle *particle_j = m_particles[j];
                if (i < j){
                    RelativeDistance = particle_i->position - particle_j->position;
                }
            }
        }
        double R12 = RelativeDistance.length();
        double R12beta = (1.0 + R12*beta);
        double R12beta2 = R12beta*R12beta;
        double LocalEnergyWithoutCoulomb = -0.5*(m_alphaomega2*R2 - 4.0*m_alphaomega - 2.0*m_a*m_alphaomega*R12/R12beta2 + (2.0*m_a/R12beta2)*(m_a/R12beta2 + 1.0/R12 - 2.0*beta/R12beta)) + 0.5*m_homega2*R2;
        if (m_InsideSteepestDescent == true){
            m_LocalEnergyWFDerivativeAlpha = -0.5*homega*R2;
            m_LocalEnergyWFDerivativeBeta = -m_a*R12*R12/R12beta2;
        }
        if (m_Coulomb == 0){
            return LocalEnergyWithoutCoulomb;
        } else {
            return LocalEnergyWithoutCoulomb + 1.0/R12;
        }
    }
}

double QuantumDot::calculateLocalEnergyWithoutJastrow(){
    double R2 = 0;
    vec3 RelativeDistance;
    if (m_Coulomb == 0) {
        for(Particle *particle : m_particles){
            double length2 = particle->position.lengthSquared();
            R2 += length2;
        }
        return 0.5*homega*homega*(1.0 - alpha*alpha)*R2 +2.0*alpha*homega;
    } else {
        for(size_t i=0; i<m_particles.size(); i++) {
            Particle *particle_i = m_particles[i];
            double length2 = particle_i->position.lengthSquared();
            R2 += length2;
            for(size_t j=0; j<m_particles.size(); j++) {
                Particle *particle_j = m_particles[j];
                if (i < j){
                    RelativeDistance = particle_i->position - particle_j->position;
                }
            }
        }
    return 0.5*m_homega2*(1.0 - m_alpha2)*R2 +2.0*m_alphaomega + 1.0/RelativeDistance.length();
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

void QuantumDot::setJastrowFactor(int trigger){
    m_Jastrow = trigger;
}

void QuantumDot::setSpinParameter(int trigger){
    m_a = trigger;
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

void QuantumDot::resetSteepestDescentHelpVars(){
    m_ExpectationLocalEnergyDerivativeAlphaSecondTerm = 0.0;
    m_ExpectationLocalEnergyDerivativeBetaSecondTerm = 0.0;
    m_ExpectationLocalEnergyDerivativeAlphaFirstTerm = 0.0;
    m_ExpectationLocalEnergyDerivativeBetaFirstTerm = 0.0;
    m_LocalEnergyWFDerivativeAlpha = 0.0;
    m_LocalEnergyWFDerivativeBeta = 0.0;
    m_MeanLocalEnergyWFDerivativeAlpha = 0.0;
    m_MeanLocalEnergyWFDerivativeBeta = 0.0;
    m_ExpectationLocalEnergyDerivativeAlpha = 0.0;
    m_ExpectationLocalEnergyDerivativeBeta = 0.0;
}


void QuantumDot::applySteepestDescent(int MonteCarloSamplesVariational,
                                      int MaxSteepestDescentIterations,
                                      double SteepestDescentStep,
                                      double tolerance){
    vec3 VarParametersOld;
    vec3 VarParametersNew;
    vec3 LocalEnergyExpectDerivative;
    m_InsideSteepestDescent = true;
    int i = 0;
    double Diff = 0.1;
    VarParametersOld.set(alpha, beta, 0.0);
    while (i < MaxSteepestDescentIterations || Diff < tolerance ){
        applyVMC(MonteCarloSamplesVariational);
        LocalEnergyExpectDerivative.set(m_ExpectationLocalEnergyDerivativeAlpha, m_ExpectationLocalEnergyDerivativeBeta, 0.0);
        VarParametersNew = VarParametersOld - LocalEnergyExpectDerivative*SteepestDescentStep;
        setVariationalParameters(VarParametersNew[0], VarParametersNew[1]);
        VarParametersOld.set(VarParametersNew[0], VarParametersNew[1],0.0);
        resetSteepestDescentHelpVars();
        i++;
        cout << "Steepest Descent iteration # " << i << endl;
        cout << "=================================" << endl;

    }
    m_InsideSteepestDescent = false;

}

void QuantumDot::setMPIenv(int argc, char **argv){
    m_argc = argc ;
    m_argv = argv;
}

