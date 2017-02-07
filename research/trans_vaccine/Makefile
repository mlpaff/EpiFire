CFLAGS=--ansi --pedantic -O2
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

all: epifire gsl vax_sim

epifire: 
	$(MAKE) -C ../../src/

gsl:
	$(MAKE) -C ../../gsl_subset/

vax_sim: vax_sim.cpp
	g++ $(CFLAGS) vax_sim.cpp $(INCLUDE) $(LDFLAGS) -o vax_sim

clean:
	rm -f vax_sim
