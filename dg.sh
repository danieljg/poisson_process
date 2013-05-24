#!/bin/bash
ifort -o dg fkdg.f90 aux.f90
./dg
