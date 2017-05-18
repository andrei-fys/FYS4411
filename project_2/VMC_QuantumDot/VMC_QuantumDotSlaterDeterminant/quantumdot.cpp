#include "quantumdot.h"
//#include "quantumforce.h"
#include <iostream>
#include <cmath>
#include <iomanip>  //on Mac setprecision
#include <fstream>
//#include "mpi.h"

using namespace std;

QuantumDot::QuantumDot(double h_omega, int ParticlesNumber){
    setUpCoordinatesCartesian(ParticlesNumber);
    setUpStatesCartesian(ParticlesNumber);
    homega = h_omega;
    m_sqrtomega = sqrt(homega);
    m_homega2 = h_omega*h_omega;
    polyRefArray[0] = &QuantumDot::H0;
    polyRefArray[1] = &QuantumDot::H1;
    polyRefArray[2] = &QuantumDot::H2;
    polyRefArray[3] = &QuantumDot::H3;
    FirstDerivRefArrayX[0] = &QuantumDot::Fi0DerivativeX;
    FirstDerivRefArrayX[1] = &QuantumDot::Fi1DerivativeX;
    FirstDerivRefArrayX[2] = &QuantumDot::Fi2DerivativeX;
    FirstDerivRefArrayX[3] = &QuantumDot::Fi3DerivativeX;
    FirstDerivRefArrayY[0] = &QuantumDot::Fi0DerivativeY;
    FirstDerivRefArrayY[1] = &QuantumDot::Fi1DerivativeY;
    FirstDerivRefArrayY[2] = &QuantumDot::Fi2DerivativeY;
    FirstDerivRefArrayY[3] = &QuantumDot::Fi3DerivativeY;
    SecondDerivRefArray[0] = &QuantumDot::Fi0Derivative2;
    SecondDerivRefArray[1] = &QuantumDot::Fi1Derivative2;
    SecondDerivRefArray[2] = &QuantumDot::Fi2Derivative2;
    SecondDerivRefArray[3] = &QuantumDot::Fi3Derivative2;
    AlphaDerivativeRefArray[0]=&QuantumDot::Fi0AlphaDerivative;
    AlphaDerivativeRefArray[1]=&QuantumDot::Fi1AlphaDerivative;
    AlphaDerivativeRefArray[2]=&QuantumDot::Fi2AlphaDerivative;
    AlphaDerivativeRefArray[3]=&QuantumDot::Fi3AlphaDerivative;
    AlphaDerivativeRefArray[4]=&QuantumDot::Fi4AlphaDerivative;
    AlphaDerivativeRefArray[5]=&QuantumDot::Fi5AlphaDerivative;
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

void QuantumDot::setUpCoordinatesCartesian(int ParticlesNumber){
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
    for(int i=0; i<ParticlesNumber; i++) {
        Particle *particle = new Particle();
        double x = RandomNumberGenerator(gen); // random number in the interval [0,1]
        double y = RandomNumberGenerator(gen);
        double z = 0;
        if (i < ParticlesNumber/2){
            particle->spin = 1;
        } else {
            particle->spin = -1;
        }
        particle->position.set(x,y,z);
        m_particles.push_back(particle);
    }
}

void QuantumDot::setUpStatesCartesian(int ParticlesNumber) {
    QuantumState m_q_state;
    int nx, ny, Ne;
        for (int i = 1; i < ParticlesNumber + 1; i++){
            for (int j = 0; j < i; j++){
                nx = j;
                ny = i - nx - 1;
                m_q_state.set(nx, ny, m_sm, m_s);
                Ne = (i + 1)*(nx + ny + 1);
                if (Ne <= ParticlesNumber) {
                    m_shells.push_back(m_q_state);
                }
            }
        }
}

void QuantumDot::setCoulombInterraction(int trigger){
    m_Coulomb = trigger;
}

void QuantumDot::setJastrowFactor(int trigger){
    m_Jastrow = trigger;
}

void QuantumDot::getQuantumDotParticlesCoordinates(){
    for(Particle *particle : m_particles){
        particle->position.print();
        cout << "Spin: " << particle->spin << endl;
    }
}

void QuantumDot::getQuantumDotStates(){
    for(QuantumState quantum_state : m_shells){
        cout << "nx = " << quantum_state.nx() << endl;
        cout << "ny = " <<quantum_state.ny() << endl;
        cout << "sm = " <<quantum_state.sm() << endl;
        cout << "s = " <<quantum_state.s() << endl;
        cout << "-----------------" << endl;
    }
}

void QuantumDot::setVariationalParameters(double first, double second){
    alpha = first;
    m_alpha2 = first*first;
    beta = second;
    m_alphaomega = alpha*homega;
    m_alphaomega2 = m_alphaomega*m_alphaomega;
}

QuantumDot::polynomialArray QuantumDot::getHetmitePolinomial(int i){
    return polyRefArray[i];

}

QuantumDot::FirstDerivativeArray QuantumDot::getGradientX(int i){
    return FirstDerivRefArrayX[i];
}

QuantumDot::FirstDerivativeArray QuantumDot::getGradientY(int i){
    return FirstDerivRefArrayY[i];
}

QuantumDot::SecondDerivativeArray QuantumDot::getLaplasian(int i){
    return SecondDerivRefArray[i];
}

QuantumDot::AlphaDerivativeArray QuantumDot::getDerivativeSDonAlpha(int i){
    return AlphaDerivativeRefArray[i];
}

void QuantumDot::setUpSlaterDeterminant(){
    QuantumDot::polynomialArray HermiteX;
    QuantumDot::polynomialArray HermiteY;
    m_SD_up_inverse.zeros(m_shells.size(),m_shells.size());
    m_SD_down_inverse.zeros(m_shells.size(),m_shells.size());
    m_SD_up.zeros(m_shells.size(),m_shells.size());
    m_SD_down.zeros(m_shells.size(),m_shells.size());
    for(size_t k=0; k< m_particles.size()/2 ; k++) {
        Particle *particle_spin_up = m_particles[k];
        Particle *particle_spin_down = m_particles[k+m_particles.size()/2];
        for(size_t l=0; l< m_shells.size() ; l++) {
            QuantumState qstate = m_shells[l];
            HermiteX = getHetmitePolinomial(qstate.nx());
            HermiteY = getHetmitePolinomial(qstate.ny());
            double expSD_up = exp(-0.5*m_alphaomega*particle_spin_up->position.lengthSquared());
            double expSD_down = exp(-0.5*m_alphaomega*particle_spin_down->position.lengthSquared());
            m_SD_up_inverse(k,l) = expSD_up*HermiteX(m_sqrtomega*particle_spin_up->position.x())*HermiteY(m_sqrtomega*particle_spin_up->position.y());
            m_SD_down_inverse(k,l) = expSD_down*HermiteX(m_sqrtomega*particle_spin_down->position.x())*HermiteY(m_sqrtomega*particle_spin_down->position.y());
        }
    }
    m_SD_up = m_SD_up_inverse;
    m_SD_down =   m_SD_down_inverse;
    m_SD_up_inverse = m_SD_up_inverse.i();
    m_SD_down_inverse = m_SD_down_inverse.i();
}

void QuantumDot::calculateQuantumForce(size_t i){
    vec3 RelativeDistance;
    double Fx = 0.0;
    double Fy = 0.0;
    Particle *moving_particle = m_particles[i];
    QuantumDot::FirstDerivativeArray GradientX;
    QuantumDot::FirstDerivativeArray GradientY;

    double SummGradX = 0.0;
    double SummGradY = 0.0;
    double exponent = exp(-0.5*m_alphaomega*moving_particle->position.lengthSquared());
    for(size_t j=0; j< m_shells.size() ; j++) {
        GradientX = getGradientX(j);
        GradientY = getGradientY(j);
        if (moving_particle->spin != 1) {
            SummGradX += m_SD_down_inverse(j,i-m_shells.size())*GradientX(moving_particle->position.x(), moving_particle->position.y(), exponent, m_alphaomega, m_sqrtomega);
            SummGradY += m_SD_down_inverse(j,i-m_shells.size())*GradientY(moving_particle->position.x(), moving_particle->position.y(), exponent, m_alphaomega, m_sqrtomega );
        } else {
            SummGradX += m_SD_up_inverse(j,i)*GradientX(moving_particle->position.x(), moving_particle->position.y(), exponent, m_alphaomega, m_sqrtomega);
            SummGradY += m_SD_up_inverse(j,i)*GradientY(moving_particle->position.x(), moving_particle->position.y(), exponent, m_alphaomega, m_sqrtomega );
        }
    }
    double abetaterm = 0.0;
    double firstsumX = 0.0;
    double firstsumY = 0.0;
    double secondsumX = 0.0;
    double secondsumY = 0.0;
    if (m_Jastrow == 1) {
        for(size_t k=0; k<i ; k++) {
            Particle *particle = m_particles[k];
            RelativeDistance = particle->position - moving_particle->position;
            if (moving_particle->spin != particle->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            }
            double reldistfirstX = (moving_particle->position.x() - particle->position.x())/RelativeDistance.length();
            double reldistfirstY = (moving_particle->position.y() - particle->position.y())/RelativeDistance.length();
            firstsumX += reldistfirstX*abetaterm;
            firstsumY += reldistfirstY*abetaterm;
        }
        for(size_t k=i+1; k< m_particles.size() ; k++) {
            Particle *particle = m_particles[k];
            RelativeDistance = particle->position - moving_particle->position;
            if (moving_particle->spin != particle->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            }
            double reldistfirstX = (particle->position.x() - moving_particle->position.x())/RelativeDistance.length();
            double reldistfirstY = (particle->position.y() - moving_particle->position.y())/RelativeDistance.length();
            secondsumX += reldistfirstX*abetaterm;
            secondsumY += reldistfirstY*abetaterm;
        }
    }
    Fx = (firstsumX - secondsumX)*2.0 +SummGradX*2.0;
    Fy = (firstsumY - secondsumY)*2.0 +SummGradY*2.0;
    m_QForceOld.setX(Fx);
    m_QForceOld.setY(Fy);
}

void QuantumDot::calculateQuantumForceNew(size_t i){
    vec3 RelativeDistance;
    double Fx = 0.0;
    double Fy = 0.0;
    Particle *moving_particle = m_particles[i];
    QuantumDot::FirstDerivativeArray GradientX;
    QuantumDot::FirstDerivativeArray GradientY;
    double SummGradX = 0.0;
    double SummGradY = 0.0;
    double exponent = exp(-0.5*m_alphaomega*moving_particle->positionNew.lengthSquared());
    for(size_t j=0; j< m_shells.size() ; j++) {
        GradientX = getGradientX(j);
        GradientY = getGradientY(j);
        if (moving_particle->spin != 1) {
            SummGradX += m_SD_down_inverse(j,i-m_shells.size())*GradientX(moving_particle->positionNew.x(), moving_particle->positionNew.y(), exponent, m_alphaomega, m_sqrtomega);
            SummGradY += m_SD_down_inverse(j,i-m_shells.size())*GradientY(moving_particle->positionNew.x(), moving_particle->positionNew.y(), exponent, m_alphaomega, m_sqrtomega );
        } else {
            SummGradX += m_SD_up_inverse(j,i)*GradientX(moving_particle->positionNew.x(), moving_particle->positionNew.y(), exponent, m_alphaomega, m_sqrtomega);
            SummGradY += m_SD_up_inverse(j,i)*GradientY(moving_particle->positionNew.x(), moving_particle->positionNew.y(), exponent, m_alphaomega, m_sqrtomega );
        }

    }
    double abetaterm = 0.0;
    double firstsumX = 0.0;
    double firstsumY = 0.0;
    double secondsumX = 0.0;
    double secondsumY = 0.0;
    if (m_Jastrow == 1) {
        for(size_t k=0; k<i ; k++) {
            Particle *particle = m_particles[k];
            RelativeDistance = particle->position - moving_particle->positionNew;
            if (moving_particle->spin != particle->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            }
            double reldistfirstX = (moving_particle->positionNew.x() - particle->position.x())/RelativeDistance.length();
            double reldistfirstY = (moving_particle->positionNew.y() - particle->position.y())/RelativeDistance.length();
            firstsumX += reldistfirstX*abetaterm;
            firstsumY += reldistfirstY*abetaterm;
        }
        for(size_t k=i+1; k< m_particles.size() ; k++) {
            Particle *particle = m_particles[k];
            RelativeDistance = particle->position - moving_particle->positionNew;
            if (moving_particle->spin != particle->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            }
            double reldistfirstX = (particle->position.x() - moving_particle->positionNew.x())/RelativeDistance.length();
            double reldistfirstY = (particle->position.y() - moving_particle->positionNew.y())/RelativeDistance.length();
            secondsumX += reldistfirstX*abetaterm;
            secondsumY += reldistfirstY*abetaterm;
        }
    }
    Fx = (firstsumX - secondsumX)*2.0  + SummGradX*2.0/(m_RSD*m_RJ);
    Fy = (firstsumY - secondsumY)*2.0  + SummGradY*2.0/(m_RSD*m_RJ);
    m_QForceNew.setX(Fx);
    m_QForceNew.setY(Fy);
}

double QuantumDot::calculateGreenFunctionRatio(size_t j){
    Particle *particle = m_particles[j];
    double GreenFun = 0.0;
    for (int i=0; i<2; i++){
        GreenFun += 0.5*(m_QForceOld[i] + m_QForceNew[i])*(D*dt*0.5*(m_QForceOld[i] - m_QForceNew[i]) - particle->positionNew[i] + particle->position[i]);
    }
    return exp(GreenFun);
}

double QuantumDot::calculateJastrowRatio(size_t i){
    if (m_Jastrow == 0) {
        m_RJ = 1.0;
        return 1.0;
    } else {
        Particle *moving_particle = m_particles[i];
        vec3 RelativeDistance;
        vec3 RelativeDistanceNew;
        double betaRold = 0.0;
        double betaRnew = 0.0;
        double aRold = 0.0;
        double aRnew = 0.0;
        double firstterm = 0.0;
        double secondterm = 0.0;
        for(size_t k=0; k<i ; k++) {
            Particle *particle = m_particles[k];
            RelativeDistanceNew = particle->position - moving_particle->positionNew;
            RelativeDistance = particle->position - moving_particle->position;
            betaRold = 1.0 + beta*RelativeDistance.length();
            betaRnew = 1.0 + beta*RelativeDistanceNew.length();
            if (moving_particle->spin != particle->spin) {
                aRold = RelativeDistance.length()*1.0;
                aRnew = RelativeDistanceNew.length()*1.0;
            } else {
                aRold = RelativeDistance.length()/3.0;
                aRnew = RelativeDistanceNew.length()/3.0;
            }
            firstterm += (aRnew/betaRnew - aRold/betaRold);
        }
        for(size_t k=i+1; k< m_particles.size() ; k++) {
            Particle *particle = m_particles[k];
            RelativeDistanceNew = particle->position - moving_particle->positionNew;
            RelativeDistance = particle->position - moving_particle->position;
            betaRold = 1.0 + beta*RelativeDistance.length();
            betaRnew = 1.0 + beta*RelativeDistanceNew.length();
            if (moving_particle->spin != particle->spin) {
                aRold = RelativeDistance.length()*1.0;
                aRnew = RelativeDistanceNew.length()*1.0;
            } else {
                aRold = RelativeDistance.length()/3.0;
                aRnew = RelativeDistanceNew.length()/3.0;
            }
            secondterm += (aRnew/betaRnew - aRold/betaRold);
        }
        double RJ = exp(firstterm + secondterm);
        m_RJ = RJ;
        return RJ;
    }
}

double QuantumDot::calculateSDRatio(size_t i){
    Particle *moving_particle = m_particles[i];
    QuantumDot::polynomialArray HermiteX;
    QuantumDot::polynomialArray HermiteY;
    double SDelement = 0.0;
    double RSD = 0.0;
    double exponent = exp(-0.5*m_alphaomega*moving_particle->positionNew.lengthSquared());
    for(size_t j=0; j< m_shells.size() ; j++) {
        QuantumState qstate = m_shells[j];
        HermiteX = getHetmitePolinomial(qstate.nx());
        HermiteY = getHetmitePolinomial(qstate.ny());
        if (moving_particle->spin == 1){
            SDelement = m_SD_up_inverse(j,i);
        } else {
            SDelement = m_SD_down_inverse(j,i-m_shells.size());
        }
        RSD += SDelement*exponent*HermiteX(m_sqrtomega*moving_particle->positionNew.x())*HermiteY(m_sqrtomega*moving_particle->positionNew.y());
    }
    m_RSD = RSD;
    return RSD;
}

void QuantumDot::updateSlaterDeterminant(size_t i){
    QuantumDot::polynomialArray HermiteX;
    QuantumDot::polynomialArray HermiteY;
    //arma::mat D_rnew(m_shells.size(),m_shells.size());
    arma::mat tempDeterminant(m_shells.size(),m_shells.size());
    Particle *moving_particle = m_particles[i];
    if (moving_particle->spin == 1){
        tempDeterminant = m_SD_up;
    } else {
        tempDeterminant = m_SD_down;
        i -= m_shells.size();
    }
    for (size_t j=0; j<m_shells.size(); j++){
        QuantumState qstate = m_shells[j];
        HermiteX = getHetmitePolinomial(qstate.nx());
        HermiteY = getHetmitePolinomial(qstate.ny());
        double expSD = exp(-0.5*m_alphaomega*moving_particle->position.lengthSquared());
        tempDeterminant(i,j) = expSD*HermiteX(m_sqrtomega*moving_particle->position.x())*HermiteY(m_sqrtomega*moving_particle->position.y());
    }

    if (moving_particle->spin == 1){
        m_SD_up = tempDeterminant;
        m_SD_up_inverse = m_SD_up.i();
    } else {
        m_SD_down = tempDeterminant;
        m_SD_down_inverse = m_SD_down.i();
    }    

}

void QuantumDot::updateInverseSlaterDeterminant(size_t i){
    QuantumDot::polynomialArray HermiteX;
    QuantumDot::polynomialArray HermiteY;
    Particle *moving_particle = m_particles[i];
    arma::mat D_inverse_rnew(m_shells.size(),m_shells.size());
    arma::mat tempDeterminant(m_shells.size(),m_shells.size());
    arma::mat tempDeterminantInv(m_shells.size(),m_shells.size());
    if (moving_particle->spin == 1){
        tempDeterminant = m_SD_up;
        tempDeterminantInv = m_SD_up_inverse;
    } else {
        tempDeterminant = m_SD_down;
        tempDeterminantInv = m_SD_down_inverse;
        i -= m_shells.size();
    }
    for(size_t k=0; k<m_shells.size(); k++) {
        for(size_t j=0; j<m_shells.size(); j++) {
            if (j == i){
                D_inverse_rnew(k,j) = tempDeterminantInv(k,i)/(m_RSD*m_RJ);
                double tempSumm = 0.0;
                for (size_t l=0; l<m_shells.size(); l++){
                    tempSumm += tempDeterminant(i,l)*tempDeterminantInv(l,j);
                }
                D_inverse_rnew(k,j) *= tempSumm;
            } else {
                D_inverse_rnew(k,j) = tempDeterminantInv(k,j);
                double tempSumm1 = 0.0;
                for (size_t l=0; l<m_shells.size(); l++){
                    QuantumState qstate = m_shells[l];
                    HermiteX = getHetmitePolinomial(qstate.nx());
                    HermiteY = getHetmitePolinomial(qstate.ny());
                    double exponent = exp(-0.5*m_alphaomega*moving_particle->positionNew.lengthSquared());
                    double DetRNew = exponent*HermiteX(m_sqrtomega*moving_particle->positionNew.x())*HermiteY(m_sqrtomega*moving_particle->positionNew.y());
                    tempSumm1 += tempDeterminantInv(l,j)*DetRNew;
                }
                tempSumm1 *= tempDeterminantInv(k,i)/(m_RSD*m_RJ);
                D_inverse_rnew(k,j) -= tempSumm1;
            }
        }
    }
    if (moving_particle->spin == 1){
        m_SD_up_inverse = D_inverse_rnew;
    } else {
        m_SD_down_inverse = D_inverse_rnew;
    }
}

double QuantumDot::calculateLocalEnergy(){
    //double LocalEnergy;

    calculateLaplasianSD(); //Main term


    if (m_Jastrow == 1 ) {  //Takes into acconunt correlations
        calculateDotProdGradientJastrowAndSD();
        calculateLaplasianJastrow();
    } else {
        m_LaplasianJastrow = 0.0;
        m_DotProdGradientJastrowAndSD = 0.0;
    }

    double HOPotentialEnergy = 0.0; // Harmonic oscillator part
    for(size_t j=0; j< m_particles.size() ; j++) {
        Particle *particle = m_particles[j];
        HOPotentialEnergy += 0.5*m_homega2*particle->position.lengthSquared();
    }

    double CoulombPotentialEnergy = 0.0; //Takes into account Coulomb interraction between electrons
    if (m_Coulomb == 1 ) {
        for(size_t i=0; i< m_particles.size() ; i++) {
            Particle *particle_i = m_particles[i];
            for(size_t j=i+1; j< m_particles.size() ; j++){
                Particle *particle_j = m_particles[j];
                CoulombPotentialEnergy += 1.0/(particle_i->position - particle_j->position).length();
            }
        }
    }
    double PotentialEnergy = 0.0;
    double KineticEnergy = 0.0;
    PotentialEnergy = HOPotentialEnergy + CoulombPotentialEnergy;
    KineticEnergy = -0.5*(m_LaplasianSD + m_LaplasianJastrow + 2.0*m_DotProdGradientJastrowAndSD);
    m_KineticEnergy += KineticEnergy;
    m_PotentialEnergy += PotentialEnergy;

    return KineticEnergy + PotentialEnergy;
}

void QuantumDot::calculateDotProdGradientJastrowAndSD(){
    m_DotProdGradientJastrowAndSD = 0.0;
    m_gradientDotProductJastrowAllParticles = 0.0;
    vec3 JastrowFactorGradient;
    vec3 RelativeDistance;
    vec3 SDGradient;
    QuantumDot::FirstDerivativeArray GradientX;
    QuantumDot::FirstDerivativeArray GradientY;
    for (size_t i=0; i<m_particles.size(); i++) {
        Particle *particle_i = m_particles[i];

        double SummGradX = 0.0;
        double SummGradY = 0.0;
        double exponent = exp(-0.5*m_alphaomega*particle_i->position.lengthSquared());
        for(size_t j=0; j< m_shells.size() ; j++) {
            GradientX = getGradientX(j);
            GradientY = getGradientY(j);
            if (particle_i->spin != 1) {
                SummGradX += m_SD_down_inverse(j,i-m_shells.size())*GradientX(particle_i->position.x(), particle_i->position.y(), exponent, m_alphaomega, m_sqrtomega);
                SummGradY += m_SD_down_inverse(j,i-m_shells.size())*GradientY(particle_i->position.x(), particle_i->position.y(), exponent, m_alphaomega, m_sqrtomega );
            } else {
                SummGradX += m_SD_up_inverse(j,i)*GradientX(particle_i->position.x(), particle_i->position.y(), exponent, m_alphaomega, m_sqrtomega);
                SummGradY += m_SD_up_inverse(j,i)*GradientY(particle_i->position.x(), particle_i->position.y(), exponent, m_alphaomega, m_sqrtomega );
            }
        }

        SDGradient.setX(SummGradX);
        SDGradient.setY(SummGradY);

        double abetaterm = 0.0;
        double firstsumX = 0.0;
        double firstsumY = 0.0;
        for(size_t k=0; k<i ; k++) {
            Particle *particle_k = m_particles[k];
            RelativeDistance = particle_k->position - particle_i->position;
            if (particle_i->spin != particle_k->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            }
            double reldistfirstX = (particle_i->position.x() - particle_k->position.x())/RelativeDistance.length();
            double reldistfirstY = (particle_i->position.y() - particle_k->position.y())/RelativeDistance.length();
            firstsumX += reldistfirstX*abetaterm;
            firstsumY += reldistfirstY*abetaterm;
        }
        abetaterm = 0.0;
        double secondsumX = 0.0;
        double secondsumY = 0.0;
        for(size_t k=i+1; k< m_particles.size() ; k++) {
            Particle *particle_k = m_particles[k];
            RelativeDistance = particle_k->position - particle_i->position;
            if (particle_i->spin != particle_k->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            }
            double reldistsecondX = (particle_k->position.x() - particle_i->position.x())/RelativeDistance.length();
            double reldistsecondY = (particle_k->position.y() - particle_i->position.y())/RelativeDistance.length();
            secondsumX += reldistsecondX*abetaterm;
            secondsumY += reldistsecondY*abetaterm;
        }
        JastrowFactorGradient.setX(firstsumX - secondsumX);
        JastrowFactorGradient.setY(firstsumY - secondsumY);
        m_DotProdGradientJastrowAndSD += (SDGradient[0]*JastrowFactorGradient[0] + SDGradient[1]*JastrowFactorGradient[1]);
    }
    m_gradientDotProductJastrowAllParticles = JastrowFactorGradient.lengthSquared();
}

void QuantumDot::calculateLaplasianJastrow(){
    vec3 RelativeDistance;
    double LaplasianJastrow = 0.0;
    for (size_t i=0; i<m_particles.size(); i++) {
        Particle *particle_i = m_particles[i];
        double abetaterm = 0.0;
        double abetaterm2 = 0.0;
        double firstsum = 0.0;
        double secondsum = 0.0;
        for(size_t k=0; k<i ; k++) {
            Particle *particle_k = m_particles[k];
            RelativeDistance = particle_k->position - particle_i->position;
            if (particle_i->spin != particle_k->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
                abetaterm2 = -2.0*beta/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
                abetaterm2 = -2.0/(3.0*((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length())));
            }
            abetaterm /= RelativeDistance.length();
            firstsum += (abetaterm + abetaterm2);
        }
        for(size_t k=i+1; k< m_particles.size() ; k++) {
            Particle *particle_k = m_particles[k];
            RelativeDistance = particle_k->position - particle_i->position;
            if (particle_i->spin != particle_k->spin) {
                abetaterm = 1.0/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
                abetaterm2 = -2.0*beta/((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
            } else {
                abetaterm = 1.0/(3.0*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length()));
                abetaterm2 = -2.0/(3.0*((1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length())*(1.0 + beta*RelativeDistance.length())));
            }
            abetaterm /= RelativeDistance.length();
            secondsum += (abetaterm + abetaterm2);
        }
    LaplasianJastrow += (firstsum + secondsum);
    }
    LaplasianJastrow += m_gradientDotProductJastrowAllParticles;
    m_LaplasianJastrow = LaplasianJastrow;
}

void QuantumDot::calculateLaplasianSD(){
    double SDLaplasianAll = 0.0;
    QuantumDot::SecondDerivativeArray Laplasian;
    for (size_t i=0; i<m_particles.size(); i++) {
        Particle *particle_i = m_particles[i];
        double SummLaplasian = 0.0;
        double exponent = exp(-0.5*m_alphaomega*particle_i->position.lengthSquared());
        for(size_t j=0; j< m_shells.size() ; j++) {
            Laplasian = getLaplasian(j);
            if (particle_i->spin != 1) {
                SummLaplasian += m_SD_down_inverse(j,i-m_shells.size())*Laplasian(particle_i->position.x(), particle_i->position.y(), exponent, m_alphaomega, m_sqrtomega);
            } else {
                SummLaplasian += m_SD_up_inverse(j,i)*Laplasian(particle_i->position.x(), particle_i->position.y(), exponent, m_alphaomega, m_sqrtomega);
            }
        }
        SDLaplasianAll += SummLaplasian;
    }
    m_LaplasianSD = SDLaplasianAll;
}

void QuantumDot::applyVMC(int MCSamples){
    setUpSlaterDeterminant();
    vector<double> LocalEnergyVector;
    //string string_rank = to_string(my_rank);
    string outputfile = "LocalEnergy_";
    //outputfile.append(string_rank);
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    normal_distribution<double> GaussianBlur(0.0, 1.0);
    m_KineticEnergy = 0.0;
    m_PotentialEnergy = 0.0;
    int MC_counter = 0;
    int accept=0;
    double MeanLocalEnergy = 0.0;
    double MeanLocalEnergy2 = 0.0;
    for(int i=0; i<MCSamples; i++){
        for(size_t j=0; j< m_particles.size() ; j++) {
            Particle *particle = m_particles[j];
            calculateQuantumForce(j);
            particle->positionNew.setX(particle->position.x() + GaussianBlur(gen)*sqrt(dt) + m_QForceOld.x()*dt*D);
            particle->positionNew.setY(particle->position.y() + GaussianBlur(gen)*sqrt(dt) + m_QForceOld.y()*dt*D);
            double SDRatio = calculateSDRatio(j);
            double JRatio = calculateJastrowRatio(j);
            calculateQuantumForceNew(j);
            double GRatio = calculateGreenFunctionRatio(j);
            double w = SDRatio*JRatio*SDRatio*JRatio*GRatio;
            double r = RandomNumberGenerator(gen);
            if (w >= r) {
               particle->position.setX(particle->positionNew.x());
               particle->position.setY(particle->positionNew.y());
               updateSlaterDeterminant(j);
               accept++;
            }
        }
        double Elocal = calculateLocalEnergy();
        double Elocal2 = Elocal*Elocal;
        MeanLocalEnergy += Elocal;
        MeanLocalEnergy2 += Elocal2;
        LocalEnergyVector.push_back(Elocal);
        if (i % 100000 == 0){ //every 10 000 000 (ten millions) MC samples
            writeVectorToBinaryFile(outputfile, LocalEnergyVector);
            LocalEnergyVector.clear();
        }
        MC_counter++;
    }
    MeanLocalEnergy /= (double) MCSamples;
    MeanLocalEnergy2 /= (double) MCSamples;
    cout << "Accept " << ((double)accept/(double)(m_particles.size()*MCSamples))*100.0 << endl;
    cout << "Elocal " << MeanLocalEnergy << endl;
    cout << "Variance " << (MeanLocalEnergy2 - MeanLocalEnergy*MeanLocalEnergy)/MCSamples << endl;
    cout << "Kinetic " << (double) m_KineticEnergy/MCSamples << endl;
    cout << "Potential " << (double) m_PotentialEnergy/MCSamples << endl;
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

void QuantumDot::steepestDescentCalculateWFderivativeOnVarParameters(){
    QuantumDot::AlphaDerivativeArray DerivativeOnAlpha;
    arma::mat MatDerivativeOnAlphaSpinUp;
    arma::mat MatDerivativeOnAlphaSpinDown;
    MatDerivativeOnAlphaSpinUp.zeros(m_shells.size(),m_shells.size());
    MatDerivativeOnAlphaSpinDown.zeros(m_shells.size(),m_shells.size());
    for(size_t k=0; k< m_particles.size()/2 ; k++) {
        Particle *particle_spin_up = m_particles[k];
        Particle *particle_spin_down = m_particles[k+m_particles.size()/2];
        for(size_t l=0; l< m_shells.size() ; l++) {
            DerivativeOnAlpha = getDerivativeSDonAlpha(l);
            double expSD_up = exp(-0.5*m_alphaomega*particle_spin_up->position.lengthSquared());
            double expSD_down = exp(-0.5*m_alphaomega*particle_spin_down->position.lengthSquared());
            MatDerivativeOnAlphaSpinUp(k,l) = DerivativeOnAlpha(particle_spin_up->position.x(), particle_spin_up->position.y(), expSD_up, m_sqrtomega, homega);
            MatDerivativeOnAlphaSpinDown(k,l) = DerivativeOnAlpha(particle_spin_down->position.x(), particle_spin_down->position.y(), expSD_down, m_sqrtomega, homega);
        }
    }
    m_LocalEnergyWFDerivativeAlpha = arma::trace(MatDerivativeOnAlphaSpinUp*m_SD_up_inverse) + arma::trace(MatDerivativeOnAlphaSpinDown*m_SD_down_inverse);
    m_LocalEnergyWFDerivativeBeta = 0.0;
    for(size_t i=0; i< m_particles.size() ; i++){
        Particle *particle_i = m_particles[i];
        for(size_t j=i+1; j< m_particles.size() ; j++){
            Particle *particle_j = m_particles[j];
            if (particle_j->spin != particle_i->spin) {
                m_LocalEnergyWFDerivativeBeta += 1.0/(1.0/(particle_i->position - particle_j->position).length() + beta)*(1.0/(particle_i->position - particle_j->position).length() + beta);
            } else {
                m_LocalEnergyWFDerivativeBeta += 1.0/(3.0*(1.0/(particle_i->position - particle_j->position).length() + beta)*(1.0/(particle_i->position - particle_j->position).length() + beta));
            }
        }
    }
}

void QuantumDot::applyVMCSteepestDescent(int MCSamples){
    setUpSlaterDeterminant();
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    normal_distribution<double> GaussianBlur(0.0, 1.0);
    m_KineticEnergy = 0.0;
    m_PotentialEnergy = 0.0;
    int MC_counter = 0;
    int accept=0;
    double ElocalSum = 0.0;
    for(int i=0; i<MCSamples; i++){
        for(size_t j=0; j< m_particles.size() ; j++) {
            Particle *particle = m_particles[j];
            calculateQuantumForce(j);
            particle->positionNew.setX(particle->position.x() + GaussianBlur(gen)*sqrt(dt) + m_QForceOld.x()*dt*D);
            particle->positionNew.setY(particle->position.y() + GaussianBlur(gen)*sqrt(dt) + m_QForceOld.y()*dt*D);
            double SDRatio = calculateSDRatio(j);
            double JRatio = calculateJastrowRatio(j);
            calculateQuantumForceNew(j);
            double GRatio = calculateGreenFunctionRatio(j);
            double w = SDRatio*JRatio*GRatio;
            double r = RandomNumberGenerator(gen);
            if (w >= r) {
               particle->position.setX(particle->positionNew.x());
               particle->position.setY(particle->positionNew.y());
               updateSlaterDeterminant(j);
               accept++;
            }
        }
        MC_counter++;
        double LocalEnergy = calculateLocalEnergy();
        steepestDescentCalculateWFderivativeOnVarParameters();
        m_MeanLocalEnergyWFDerivativeAlpha += m_LocalEnergyWFDerivativeAlpha;
        m_MeanLocalEnergyWFDerivativeBeta += m_LocalEnergyWFDerivativeBeta;
        m_ExpectationLocalEnergyDerivativeAlphaFirstTerm += m_LocalEnergyWFDerivativeAlpha*LocalEnergy;
        m_ExpectationLocalEnergyDerivativeBetaFirstTerm += m_LocalEnergyWFDerivativeBeta*LocalEnergy;
        ElocalSum += LocalEnergy;
    }
    double Energy = (double) ElocalSum/MCSamples;
    m_ExpectationLocalEnergyDerivativeAlphaSecondTerm = Energy*((double)m_MeanLocalEnergyWFDerivativeAlpha/MCSamples);
    m_ExpectationLocalEnergyDerivativeBetaSecondTerm = Energy*((double)m_MeanLocalEnergyWFDerivativeBeta/MCSamples);
    m_ExpectationLocalEnergyDerivativeAlphaFirstTerm /= (double)MCSamples;
    m_ExpectationLocalEnergyDerivativeBetaFirstTerm /= (double)MCSamples;
    m_ExpectationLocalEnergyDerivativeAlpha = 2.0*(m_ExpectationLocalEnergyDerivativeAlphaFirstTerm-m_ExpectationLocalEnergyDerivativeAlphaSecondTerm);
    m_ExpectationLocalEnergyDerivativeBeta = 2.0*(m_ExpectationLocalEnergyDerivativeBetaFirstTerm-m_ExpectationLocalEnergyDerivativeBetaSecondTerm);
    cout << "Accept " << ((double)accept/(double)(m_particles.size()*MCSamples))*100.0 << endl;
    cout << "Elocal " << Energy << endl;
}

void QuantumDot::applySteepestDescent(int MonteCarloSamplesVariational,
                                      int MaxSteepestDescentIterations,
                                      double SteepestDescentStep,
                                      double tolerance){
    vec3 VarParametersOld;
    vec3 VarParametersNew;
    vec3 LocalEnergyExpectDerivative;
    int i = 0;
    double Diff = 0.1;
    VarParametersOld.set(alpha, beta, 0.0);
    while (i < MaxSteepestDescentIterations && Diff < tolerance ){
        cout << "alpha " << alpha << endl;
        cout << "beta " << beta << endl;
        applyVMCSteepestDescent(MonteCarloSamplesVariational);
        LocalEnergyExpectDerivative.set(m_ExpectationLocalEnergyDerivativeAlpha, m_ExpectationLocalEnergyDerivativeBeta, 0.0);
        VarParametersNew = VarParametersOld - LocalEnergyExpectDerivative*SteepestDescentStep;
        setVariationalParameters(VarParametersNew[0], VarParametersNew[1]);
        VarParametersOld.set(VarParametersNew[0], VarParametersNew[1],0.0);
        resetSteepestDescentHelpVars();
        i++;
        cout << "Steepest Descent iteration # " << i << endl;
        cout << "=================================" << endl;

    }
}


void QuantumDot::writeVectorToBinaryFile(string ResultsFile, vector<double>& Vector){
    ofstream ofile;
    ofile.open(ResultsFile, ios::app | ios::binary);
    ofile.write(reinterpret_cast<const char*>(&Vector[0]), Vector.size() * sizeof(double));
    //ofile << setprecision(12);
    //for(double element : Vector){
    //    ofile. << element << endl;
    //}
    ofile.close();
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

/*
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

*/

/*



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

*/
