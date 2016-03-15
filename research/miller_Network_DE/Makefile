CFLAGS=--ansi --pedantic -O2
#CFLAGS=--ansi --pedantic -g
INCLUDE= -I../src/
LDFLAGS=  ../src/*.o

all: epifire gsl miller_2013_fig3 nld_fig_run try_networks

epifire:
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/

miller_2013_fig3: miller_2013_fig3.cpp
	g++ $(CFLAGS) miller_2013_fig3.cpp $(INCLUDE) -I../gsl_subset/ $(LDFLAGS) ../gsl_subset/*.o -o miller_2013_fig3

nld_fig_run: nld_fig_run.cpp
	g++ $(CFLAGS) nld_fig_run.cpp $(INCLUDE) -I../gsl_subset/ $(LDFLAGS) ../gsl_subset/*.o -o nld_fig_run

try_networks: try_networks.cpp
	g++ $(CFLAGS) try_networks.cpp $(INCLUDE) -I../gsl_subset/ $(LDFLAGS) ../gsl_subset/*.o -o try_networks

clean:
	rm -f miller_2013_fig3 nld_fig_run try_networks


