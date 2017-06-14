How to compile the source

#gfortran -O3 -o fw.exe fw.f90 -llapack -I/usr/include -lfftw3 

Execute the compiled code and produce wafefunction data 

#./fw.exe > plot.dat

Plot the data (press any key in terminal windows in order to quit the gnuplot)

#gnuplot fwplot.gnu 

Plot will be saved into PNG file as plot.png



