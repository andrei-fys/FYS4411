#!/bin/bash

#set -x

progname="Hartree-Fock_QuantumDot"
omega=1
NumberOfShells=(3 4 5 6 7 8 9 10 11 12 13 14)

function start_for_different_shells {
    for NumOfShells in ${NumberOfShells[@]}
        do
           ./"$progname" $NumOfShells $1 $omega
        done
}


NumberOfElectrons=(2 6 12 20)

for i in ${NumberOfElectrons[@]}
    do
       start_for_different_shells $i
    done

