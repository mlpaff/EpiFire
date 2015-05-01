CFLAGS=--ansi --pedantic -O2
#CFLAGS=--ansi --pedantic -g
INCLUDE= -I../src/
LDFLAGS=  ../src/*.o

all: epifire gsl miller_2013_fig3

epifire:
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/

miller_2013_fig3: miller_2013_fig3.cpp
	g++ $(CFLAGS) miller_2013_fig3.cpp $(INCLUDE) -I../gsl_subset/ $(LDFLAGS) ../gsl_subset/*.o -o miller_2013_fig3

clean:
	rm -f miller_2013_fig3


