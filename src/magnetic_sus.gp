set term jpeg
set output "emag_sus.jpg"
set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind
set object 1 rect fc  rgb "white"  fillstyle solid 1.0

set autoscale
set nokey
set grid
set title "Magnetic Susceptability of Ising Model"
set xlabel "beta"
set ylabel "X"
plot "example2.dat" using 1:6:7 with errorbars, "example2.dat" using 1:6 with lines
