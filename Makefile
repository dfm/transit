CFLAGS = -Iinclude -I/usr/local/include

default: ldtest

.cpp.o:
	g++ ${CFLAGS} -o $*.o -c $*.cpp

ldtest: src/ld_test.o src/quad.o
	g++ ${CFLAGS} src/ld_test.o src/quad.o -o bin/ldtest

clean:
	rm -rf src/*.o bin/*
