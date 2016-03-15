CFLAGS=--ansi --pedantic -O2 --std=c++11
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../src/
LDFLAGS=  ../src/*.o

default: two_percolation headstart_percolation

all: epifire gsl two_percolation

epifire: 
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/

two_percolation: two_percolation.cpp ../src/Two_Strain_Percolation_Sim.h
	g++ $(CFLAGS) two_percolation.cpp $(INCLUDE) $(LDFLAGS) -o two_percolation

headstart_percolation: headstart_percolation.cpp ../src/Two_Strain_Percolation_Sim.h
	g++ $(CFLAGS) headstart_percolation.cpp $(INCLUDE) $(LDFLAGS) -o headstart_percolation

clean:
	rm -f two_percolation
