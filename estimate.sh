#!/bin/bash
echo ================================================
echo Driver for the estimation of phase and constrast
echo and the characterization of their errors
echo By Daniel Jimenez
echo ================================================
echo Compiling
source compile.sh
echo ================================================
echo Setting up variables and cleaning previous files
rm results/phase_error*.dat
# number of fringes for each configuration
number_of_fringes=1000
# number of expected photon arrivals goes
# from : 10*2^jmin (should be 1)
# to   : 10*2^jmax (should be 10 for final test)
# (calculation is internal)
 jmin=1
 jmax=10
# evaluated visibilities go as
# vis=k*5 (internal division by 100)
# (should be from 1 to 20)
 kmin=1
 kmax=19
echo jmin=$jmin, jmax=$jmax
echo kmin=$kmin, kmax=$kmax
echo ================================================
echo main loop
# set expected number of photons
for (( j=jmin ; j<=jmax ; j++ ))
do
 tar -zxvf data/nbar_$j.tar.gz -C data/ > /dev/null
 echo -----------------------------------------------
 echo photon number loop: j=$j of $jmax
 mkdir results/nbar_$j/
 # set visibility
 for (( k=kmin ; k<=kmax ; k++))
 do
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  echo contrast loop: k=$k of $kmax
  vis=$((5*$k))
  mkdir results/nbar_$j/vis_$vis/
  ./estimate $j $vis >> results/phase_error_$vis.dat
  mv phase_histogram.dat results/nbar_$j/vis_$vis/
  mv counts_histogram.dat results/nbar_$j/vis_$vis/
  mv phase_moments.dat results/nbar_$j/vis_$vis/
  #mv fact_moments.dat results/nbar_$j/vis_$vis/
 done
 rm -r data/nbar_$j/
done

echo experimento terminado
play -n synth 1 sine 500-400
