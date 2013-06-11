#!/bin/bash
cd source
ifort -c aux.f90 estimate_module.f90 -i8 -r8 -fast -xHost -parallel
ifort -o estimate estimate_module.f90 aux.f90 mkl_dfti.f90 estimate.f90 -mkl -i8 -r8 -fast -xHost -parallel
ifort -c aux.f90 -mkl -fast -xHost
ifort -o dg_cli aux.f90 mkl_dfti.f90 dg_cli.f90 -mkl -fast -xHost
mv estimate ../
mv dg_cli ../
cd ..
