#!/usr/bin/gnuplot

reset

set terminal pdf transparent
set output 'fg-performance-feval-1-3.pdf'

# png
# set terminal pngcairo transparent size 800,400 enhanced font 'Verdana,9'
# set output 'farfield-memory.png'
# svg
#set terminal svg size 410,250 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '9' rounded dashed
#set output 'nice_web_plot.svg'

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror
set xtics (16384, 32768, 65536, 131072)

# define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

set xrange [16384:131072]
# color definitions
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 3 lc rgb '#800080' pt 8 ps 1 lt 1 lw 2 # --- purple
set style line 4 lc rgb '#00479D' pt 10 ps 1 lt 1 lw 2 # --- blue

# place legend
set key top left

# set labels
set xlabel 'Dimension der Matrix'
set ylabel 'Rechenzeit in ms'

set datafile separator ","

plot 'performance-feval-1-3.csv' u 1:2 t '1. Phase' w lp ls 1, \
     ''                          u 1:3 t '2. Phase' w lp ls 3, \
     ''                          u 1:4 t '3. Phase' w lp ls 2

#plot 'memory-galerkin-h2.csv' u 1:2 t 'Galerkin'  w ls 1, \
#     ''                       u 1:3 t 'H2-Matrix' w ls 2
