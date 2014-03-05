.cpp.o:
	g++ -Iinclude -I/usr/local/include/eigen3 -o $*.o -c $*.cpp

ldtest: src/ld_test.o src/quad.o
	g++ src/ld_test.o src/quad.o -o bin/ldtest

clean:
	rm -rf src/ld_test.o bin/ldtest
