#!/bin/bash
echo ================================================
echo Driver for the generation of poissonian fringes
echo By Daniel Jimenez
echo ================================================
echo Compiling
source compile.sh
echo ================================================
echo Setting up variables and cleaning previous files
# old files be-gone
rm fakedata.*
rm data/*.gz
# simple bash driver for dg_cli
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
echo $number_of fringes repetitions to be performed
echo ================================================
echo main loop
# set expected number of photons
for (( j=jmin ; j<=jmax ; j++ ))
do
 echo -----------------------------------------------
 echo photon number loop: j=$j of $jmax
 mkdir data/nbar_$j/
 # set visibility
 for (( k=kmin ; k<=kmax ; k++))
 do
  echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  echo contrast loop: k=$k of $kmax
  vis=$((5*$k))
  mkdir data/nbar_$j/vis_$vis/
  for ((l=1 ; l<=number_of_fringes ; l++))
  do
   ./dg_cli $j $vis
   gzip fakedata.dat
   mv fakedata.dat.gz data/nbar_$j/vis_$vis/$l.gz
  done
 done
 tar -zcvf data/nbar_$j.tar.gz data/nbar_$j > /dev/null
 rm -r data/nbar_$j
done
