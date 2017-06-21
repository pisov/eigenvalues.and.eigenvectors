set size square
set xrange [-5:5]
set yrange [-1:11]
set xtics 1
set xlabel "x"
set ylabel "Energy"
plot "./plot.dat" u 1:2 w l lw 1.5 notitle
replot "./plot.dat" u 1:3 w l lw 1.5 notitle
replot "./plot.dat" u 1:4 w l lw 1.5 notitle
replot "./plot.dat" u 1:6 w l lw 1.5 notitle
pause -1
set term png
set output "plot.png"
replot
