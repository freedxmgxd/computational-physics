set terminal qt
set title 'Harmonic Oscillator Wavefunctions'
set xlabel 'x'
set ylabel 'Wavefunction'
set grid
plot 'results_b.dat' using 1:2 with lines title 'n=30'
