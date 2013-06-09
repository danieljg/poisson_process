set term postscript eps enhanced color dashed

binc(x,s) = s*(int(x/s)+0.5)
set ylabel 'relative frequency'
set xlabel 'estimated phase'
set output 'histogram_phase_error_sql.eps'
set title 'Histogram of phase estimation error \n 1e5 reps, 4e3 photons/fringe, phase=512'
plot './phase.dat' u (binc($2,1.)):(1./1024.) smooth frequency with boxes ti 'Estimated phase'
set output
