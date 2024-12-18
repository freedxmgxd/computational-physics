set terminal png size 800,600
set output 'plot.png'
set title 'Linear Congruential Generator Results'
set xlabel 'Iteration'
set ylabel 'Random Value'
plot 'data.dat' using 1:2 with points title 'Random Numbers'
