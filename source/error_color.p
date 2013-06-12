set encoding utf8
set term postscript eps enhanced color size 5in,2.8in "FreeSans-Bold" 14
set output 'error_color.eps'

#set logscale x
set logscale y
set xrange [0:100]
set yrange [8:12000]

set xlabel "Contrast (%)"
set ylabel "Mean number of photons"
set title "Color map for the expected error in phase estimation"

set pm3d
set view map
unset surface

splot './phase_error.dat' u ($1*100):($2):($3*2*pi/1024) noti
set output
