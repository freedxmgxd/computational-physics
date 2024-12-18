set terminal png
set output 'energy_plot.png'
set xlabel 'Step'
set ylabel 'Energy'
plot 'eplot.dat' with lines title 'Energy'
