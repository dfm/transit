CFLAGS = -Isrc
CLIBS  = -L/usr/local/Cellar/gfortran/4.8.2/gfortran/lib -lgfortran
MA =  mandel-agol/occultquad.o

.f.o:
	gfortran -o $*.o -c $*.f

.cpp.o:
	g++ $(CFLAGS) -o $*.o -c $*.cpp

test: src/test.o src/quad.o $(MA)
	g++ $(MA) src/test.o src/quad.o $(CLIBS) -o quad
