set title 'Simulated Annealing Trajectory with Function'
set xlabel 'x'
set ylabel 'f(x)'
plot 'function_data.dat' with lines title 'Function', \
     'annealing_data.dat' with points pointtype 7 title 'Annealing'
