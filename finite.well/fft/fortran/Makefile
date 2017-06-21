FC=gfortran
CFLAGS=-O3 -I/usr/include
LDFLAGS=-llapack -lfftw3
EXE=fw.x
SRC=fw.f90
OBJ=$(SRC:.f90=.o)

all: $(SRC) $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) -o $(EXE)


%.o: %.f90 
	$(FC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.mod $(OBJ) $(EXE)
