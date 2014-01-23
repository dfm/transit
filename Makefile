CFLAGS = -Isrc

.cpp.o:
	g++ $(CFLAGS) -o $*.o -c $*.cpp

test: src/test.o src/quad.o
	g++ src/test.o src/quad.o $(CLIBS) -o quad
