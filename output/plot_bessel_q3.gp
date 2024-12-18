set logscale xy
set xlabel 'N'
set ylabel 'D_S(N)'
set title 'D_S vs N'
set grid
set key top left
plot 'output/bessel_data_q3.dat' using 1:2 with points pointtype 7 title 'D_S_1'
