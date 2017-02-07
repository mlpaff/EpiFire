CFLAGS=--ansi --pedantic -O2
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

all: epifire gsl er_vax_sim

epifire: 
	$(MAKE) -C ../../src/

gsl:
	$(MAKE) -C ../../gsl_subset/

er_vax_sim: er_vax_sim.cpp
	g++ $(CFLAGS) er_vax_sim.cpp $(INCLUDE) $(LDFLAGS) -o er_vax_sim

clean:
	rm -f er_vax_sim
