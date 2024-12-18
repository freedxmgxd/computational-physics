set logscale xy
set xlabel 'N'
set ylabel 'D_S(N)'
set title 'D_S_1 and D_S_2 vs N'
set grid
set key top left
plot 'output/bessel_data.dat' using 1:2 with points pointtype 7 title 'D_S_1', 'output/bessel_data.dat' using 1:3 with points pointtype 5 title 'D_S_2'
