#!/usr/bin/gnuplot

reset

# png
set terminal pngcairo transparent size 800,400 enhanced font 'Verdana,9'
set output 'h2-memory.png'
# svg
#set terminal svg size 410,250 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '9' rounded dashed
#set output 'nice_web_plot.svg'

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror

# define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

# color definitions
set style line 1 lc rgb "red" # --- red
set style line 2 lc rgb '#5e9c36' # --- green
set style line 3 lc rgb '#800080' # --- purple
set style line 4 lc rgb '#00479D' # --- blue

# place legend
set key top left

# set labels
set xlabel 'Dimension of the (square) matrix'
set ylabel 'Memory in GB'

set datafile separator ","

set style data histogram
set style histogram cluster gap 1
set style fill solid

set boxwidth 0.75

plot 'memory.csv' u 2:xtic(1) t 'H2-matrix' ls 1, \
     ''           u 4 t 'nearfield matrices' ls 3, \
     ''           u 3 t 'farfield coupling matrices'  ls 4
