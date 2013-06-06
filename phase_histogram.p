set term postscript eps enhanced color dashed

bin(x,s) = s*(int(x/s))
binc(x,s) = s*(int(x/s)+0.5)
set ylabel 'relative frequency'
set xlabel 'estimated phase'
set output 'fourier_mopm.eps'
set title 'Histogram of phase estimation error for 1200 photons/fringe and Vis=0.75'
plot './phase.dat' u (binc($2,1.)):(1./1024.) smooth frequency with boxes ti 'Estimated phase'

