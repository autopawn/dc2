SOURCES = src/main.c \
    ./src/bnb.c \
    ./src/construction.c \
    ./src/expand.c \
    ./src/load.c \
    ./src/localsearch.c \
    ./src/localsearch_resende.c \
    ./src/localsearch_whitaker.c \
    ./src/output.c \
    ./src/problem.c \
    ./src/redstrategy.c \
    ./src/reduction.c \
    ./src/reduction_diversity.c \
    ./src/reduction_rank.c \
    ./src/reduction_vr.c \
    ./src/rundata.c \
    ./src/runinfo.c \
    ./src/runprecomp.c \
    ./src/shuffle.c \
    ./src/solution.c \
    ./src/utils.c


SOURCES_OPT_CHECKER = src/main_opt_checker.c \
    src/redstrategy.c \
    src/utils.c \
    src/load.c \
    src/problem.c \
    src/solution.c


compile:
	rm -rf bin || true
	mkdir bin
	gcc -g -O4 -march=native -flto -Wall $(SOURCES) -lpthread -lm -o bin/dc
	gcc -g -O2 -Wall $(SOURCES) -lpthread -lm -o bin/dc_O2
	gcc -g -pg -O4 -march=native -flto -Wall $(SOURCES) -lpthread -lm -o bin/dc_prof
	gcc -g -pedantic -Wall $(SOURCES) -lpthread -lm -D DEBUG -o bin/dc_debug
	gcc -g -O4 -march=native -flto -Wall $(SOURCES_OPT_CHECKER) -lpthread -lm -o bin/opt_checker

