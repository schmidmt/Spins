set term jpeg
set output "energy.jpg"
set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind
set object 1 rect fc  rgb "white"  fillstyle solid 1.0

set autoscale
set nokey
set grid
set title "Energy of Ising Model"
set xlabel "beta"
set ylabel "<E>"
plot "example2.dat" using 1:4:5 with errorbars, "example2.dat" using 1:4 with lines
