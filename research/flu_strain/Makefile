CFLAGS=--ansi --pedantic -O2
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

all: epifire gsl gillespie_network_two_strain gillespie_plaw_two_strain gillespie_uniform gillespie_er test

epifire: 
	$(MAKE) -C ../../src/

gsl:
	$(MAKE) -C ../../gsl_subset/


gillespie_network_two_strain: gillespie_network_two_strain.cpp
	g++ $(CFLAGS) gillespie_network_two_strain.cpp $(INCLUDE) $(LDFLAGS) -o gillespie_network_two_strain

gillespie_plaw_two_strain: gillespie_plaw_two_strain.cpp
	g++ $(CFLAGS) gillespie_plaw_two_strain.cpp $(INCLUDE) $(LDFLAGS) -o gillespie_plaw_two_strain

gillespie_uniform: gillespie_uniform.cpp
	g++ $(CFLAGS) gillespie_uniform.cpp $(INCLUDE) $(LDFLAGS) -o gillespie_uniform

gillespie_er: gillespie_er.cpp
	g++ $(CFLAGS) gillespie_er.cpp $(INCLUDE) $(LDFLAGS) -o gillespie_er

test: test.cpp
	g++ $(CFLAGS) test.cpp $(INCLUDE) $(LDFLAGS) -o test

clean:
	rm -f gillespie_network_two_strain gillespie_plaw_two_strain gillespie_uniform gillespie_er test
