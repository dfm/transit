CFLAGS = -Isrc
CLIBS  = -L/usr/local/Cellar/gfortran/4.8.2/gfortran/lib -lgfortran

.f.o:
	gfortran -o $*.o -c $*.f

.cpp.o:
	g++ $(CFLAGS) -o $*.o -c $*.cpp

test: src/test.o src/quad.o src/occultquad.o
	g++ src/occultquad.o src/test.o src/quad.o $(CLIBS) -o quad
