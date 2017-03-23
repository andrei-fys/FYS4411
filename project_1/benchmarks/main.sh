#!/bin/bash


#copy this file to compiled directory and run as
# $ ./main.sh



progname="Hartree-Fock_QuantumDot"
omega=0.1
NumberOfShells=(4 5 6)

function start_for_different_shells {
    for NumOfShells in ${NumberOfShells[@]}
        do
           ./"$progname" $NumOfShells $1 $omega
           date
        done
}


NumberOfElectrons=(2 6 12 20)

for i in ${NumberOfElectrons[@]}
    do
       start_for_different_shells $i
    done

