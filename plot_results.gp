set terminal qt
set title 'Harmonic Oscillator Wavefunctions'
set xlabel 'x'
set ylabel 'Wavefunction'
set grid
plot 'results.dat' using 1:2 with lines title 'n=0', \
     'results.dat' using 1:3 with lines title 'n=1', \
     'results.dat' using 1:4 with lines title 'n=2', \
     'results.dat' using 1:5 with lines title 'n=3'
