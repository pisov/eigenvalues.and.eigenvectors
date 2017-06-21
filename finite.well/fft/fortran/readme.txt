1. How to compile the source

gfortran -O3 -o fw.x fw.f90 -llapack -I/usr/include -lfftw3 

or

make

2. Execute the compiled code and produce wafefunction data 

./fw.x > plot.dat

3. Plot the data (press any key in terminal windows in order to quit the gnuplot)

gnuplot plot.gnu 

[Remark:] Plot will be saved into PNG file as plot.png
