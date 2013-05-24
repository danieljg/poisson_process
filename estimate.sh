#!/bin/bash
clear
echo compilando
ifort -c aux.f90 -i8 -r8 -fast -xHost -parallel
ifort -o estimate aux.o estimate.f90 -i8 -r8 -fast -xHost -parallel -ipo
echo corriendo
time ./estimate
#play -n synth 0.5 sine 300-3300
