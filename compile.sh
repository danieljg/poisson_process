#!/bin/bash
ifort -o estimate estimate_module.f90 aux.f90 mkl_dfti.f90 estimate.f90 -mkl
ifort -o dg_cli aux.f90 mkl_dfti.f90 dg_cli.f90 -mkl