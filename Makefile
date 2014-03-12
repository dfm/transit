default: ldtest integrate_test

.cpp.o:
	g++ -Iinclude -o $*.o -c $*.cpp

ldtest: src/ld_test.o src/quad.o
	g++ src/ld_test.o src/quad.o -o bin/ldtest

integrate_test: src/integrate_test.o src/quad.o
	g++ src/integrate_test.o src/quad.o -o bin/integrate_test

clean:
	rm -rf src/*.o bin/*
