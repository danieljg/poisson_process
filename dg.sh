#!/bin/bash
source compile.sh
# simple bash driver for dg_cli
# number of expected photon arrivals goes
# from : 10*2^jmin (should be 1)
# to   : 10*2^jmax (should be 10 for final test)
# (calculation is internal)
 jmin=1
 jmax=4
# evaluated visibilities go as
# vis=k*5 (internal division by 100)
 kmin=4
 kmax=19
# set expected number of photons
for (( j=jmin ; j<=jmax ; j++ ))
do
 # set visibility
 for (( k=kmin; k<=kmax ; k++))
 do
  vis=$((5*$k))
  ./dg_cli $j $vis
 done
done
