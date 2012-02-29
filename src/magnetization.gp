set term svg
set output "magnetization.svg"
set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind
set object 1 rect fc  rgb "white"  fillstyle solid 1.0

set autoscale
set nokey
set grid
set title "Magnetization of 2D Ising Model"
set xlabel "T"
set ylabel "<m>"
plot "example2.dat" using 1:2:3 with errorbars, "example2.dat" using 1:2
