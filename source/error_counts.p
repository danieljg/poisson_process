set encoding utf8
set term postscript eps enhanced color dashed size 5in,2.1in "FreeSans" 22
set output 'error_counts.eps'

set style line 1 lt rgb "#A00000" lw 4 pt 1
set style line 2 lt rgb "#008000" lw 4 pt 6
set style line 3 lt rgb "#0000A0" lw 4 pt 3
set style line 4 lt rgb "#208060" lw 4 pt 4
set style line 5 lt rgb "#606000" lw 4 pt 5

set ylabel 'Expected error [degs]'
set xlabel 'Number of photons per fringe'
#set title 'Expected error for the Fourier Method'

set ytics 30
set mytics 2
set logscale x
#set logscale y
set xrange [8:12000]
set yrange [0:110]
set tmargin 0.42

plot './phase_error_5.dat' u ($2):($3*2*180/1024) w l ls 1 ti "Contrast: 0.05",\
    './phase_error_25.dat' u ($2):($3*2*180/1024) w l ls 2 ti "0.25",\
    './phase_error_50.dat' u ($2):($3*2*180/1024) w l ls 3 ti "0.50",\
    './phase_error_75.dat' u ($2):($3*2*180/1024) w l ls 4 ti "0.75",\
    './phase_error_95.dat' u ($2):($3*2*180/1024) w l ls 5 ti "0.95"
set output
