CFLAGS = -Isrc

.cpp.o:
	g++ $(CFLAGS) -o $*.o -c $*.cpp

test: src/test.o src/quad.o src/driver.o
	g++ src/test.o src/quad.o src/driver.o $(CLIBS) -o quad
