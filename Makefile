CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I NFLlib/include/ -I NFLlib/include/nfl -I NFLlib/include/nfl/prng -I include -DNFL_OPTIMIZED=ON -DNTT_AVX2
INCLUDES = include/bench.h include/cpucycles.h
BLAKE3 = src/blake3/blake3.c src/blake3/blake3_dispatch.c src/blake3/blake3_portable.c src/blake3/blake3_sse2_x86-64_unix.S src/blake3/blake3_sse41_x86-64_unix.S src/blake3/blake3_avx2_x86-64_unix.S src/blake3/blake3_avx512_x86-64_unix.S
BENCH = src/bench.c src/cpucycles.c
TEST = src/test.c
LIBS = deps/libnfllib_static.a -lgmp -lmpfr -L deps/ -lflint -lquadmath

all: bdlop bgv shuffle piaex pianex

bdlop: src/bdlop.cpp src/bgv.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/bgv.cpp -o bgv.o
	${CPP} ${CFLAGS} -DMAIN src/bdlop.cpp bgv.o ${TEST} ${BENCH} -o bdlop ${LIBS}

bgv: src/bgv.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN src/bgv.cpp ${TEST} ${BENCH} -o bgv ${LIBS}

shuffle: src/shuffle.cpp src/bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/bdlop.cpp -o bdlop.o
	${CPP} ${CFLAGS} -DMAIN src/shuffle.cpp sample_z_small.o bdlop.o ${TEST} ${BENCH} ${BLAKE3} -o shuffle ${LIBS}

piaex: src/bdlop.cpp src/piaex.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIZE=3 -c src/bdlop.cpp -o bdlop.o
	${CPP} ${CFLAGS} -DSIZE=3 -DMAIN src/piaex.cpp bdlop.o ${TEST} ${BENCH} ${BLAKE3} -o piaex ${LIBS}

pianex: src/pianex.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -DMAIN src/pianex.cpp sample_z_small.o ${TEST} ${BENCH} ${BLAKE3} -o pianex ${LIBS}

clean:
	rm *.o bdlop bgv shuffle piaex pianex
