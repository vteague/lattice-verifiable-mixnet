CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I include/nfl/prng -I include -DNFL_OPTIMIZED=ON -DNTT_AVX2
INCLUDES = include/bench.h include/cpucycles.h
BLAKE3 = src/blake3/blake3.c src/blake3/blake3_dispatch.c src/blake3/blake3_portable.c src/blake3/blake3_sse2_x86-64_unix.S src/blake3/blake3_sse41_x86-64_unix.S src/blake3/blake3_avx2_x86-64_unix.S src/blake3/blake3_avx512_x86-64_unix.S
BENCH = src/bench.c src/cpucycles.c
TEST = src/test.c
LIBS = deps/libnfllib_static.a -lgmp -lmpfr -L deps/ -lflint

all: commit encrypt shuffle piaex pianex

commit: src/commit.cpp src/encrypt.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/encrypt.cpp -o encrypt.o
	${CPP} ${CFLAGS} -DMAIN src/commit.cpp encrypt.o ${TEST} ${BENCH} -o commit ${LIBS}

encrypt: src/encrypt.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN src/encrypt.cpp ${TEST} ${BENCH} -o encrypt ${LIBS}

shuffle: src/shuffle.cpp src/commit.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c src/gaussian_ct.cpp -o gaussian.o
	${CPP} ${CFLAGS} -c src/commit.cpp -o commit.o
	${CPP} ${CFLAGS} -DMAIN src/shuffle.cpp gaussian.o commit.o ${TEST} ${BENCH} ${BLAKE3} -o shuffle ${LIBS}

piaex: src/commit.cpp src/piaex.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIZE=3 -c src/commit.cpp -o commit.o
	${CPP} ${CFLAGS} -DSIZE=3 -DMAIN src/piaex.cpp commit.o ${TEST} ${BENCH} ${BLAKE3} -o piaex ${LIBS}

pianex: src/pianex.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_ANEX -c src/gaussian_ct.cpp -o gaussian.o
	${CPP} ${CFLAGS} -DMAIN src/pianex.cpp gaussian.o ${TEST} ${BENCH} ${BLAKE3} -o pianex ${LIBS}

clean:
	rm *.o commit encrypt shuffle piaex pianex
