set terminal postscript enhanced eps color
set output "specific_heat.eps"
set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind
set object 1 rect fc  rgb "white"  fillstyle solid 1.0

set autoscale
set nokey
unset grid
set title "Specific Heat of Ising Model"
set xlabel "{/Symbol b}"
set ylabel "c"
plot "16x16.dat" using 1:8:9 with errorbars, "16x16.dat" using 1:8 with lines
