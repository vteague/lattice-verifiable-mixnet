CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I include/nfl/prng -I include -DNFL_OPTIMIZED=ON -DNTT_AVX2
INCLUDES = include/bench.h include/cpucycles.h
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
	${CPP} ${CFLAGS} -DMAIN src/shuffle.cpp src/sha224-256.c gaussian.o commit.o ${TEST} ${BENCH} -o shuffle ${LIBS}

piaex: src/commit.cpp src/piaex.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIZE=3 -c src/commit.cpp -o commit.o
	${CPP} ${CFLAGS} -DSIZE=3 -DMAIN src/piaex.cpp src/sha224-256.c commit.o ${TEST} ${BENCH} -o piaex ${LIBS}

pianex: src/pianex.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_ANEX -c src/gaussian_ct.cpp -o gaussian.o
	${CPP} ${CFLAGS} -DMAIN src/pianex.cpp src/sha224-256.c gaussian.o ${TEST} ${BENCH} -o pianex ${LIBS}

clean:
	rm *.o commit encrypt shuffle piaex pianex
