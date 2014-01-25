CFLAGS = -Isrc

.cpp.o:
	g++ $(CFLAGS) -o $*.o -c $*.cpp

test: transit/test.o transit/quad.o transit/driver.o
	g++ transit/test.o transit/quad.o transit/driver.o $(CLIBS) -o bin/test
