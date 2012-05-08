set term jpeg
set output "magnetization.jpg"
set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind
set object 1 rect fc  rgb "white"  fillstyle solid 1.0

set autoscale
set nokey
set grid
set title "Magnetization of Ising Model"
set xlabel "beta"
set ylabel "<m>"
plot "32x32.dat" using 1:2:3 with errorbars, "32x32.dat" using 1:2 with lines
