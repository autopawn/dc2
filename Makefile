SOURCES = src/main.c \
    src/bnb.c \
    src/expand.c \
    src/localsearch.c \
    src/output.c \
    src/redstrategy.c \
    src/reduction_diversity.c \
    src/reduction_vr.c \
    src/utils.c \
    src/construction.c \
    src/load.c \
    src/problem.c \
    src/reduction.c \
    src/reduction_rank.c \
    src/shuffle.c \
    src/solution.c

SOURCES_OPT_CHECKER = src/main_opt_checker.c \
    src/redstrategy.c \
    src/utils.c \
    src/load.c \
    src/problem.c \
    src/solution.c


compile:
	rm -rf bin || true
	mkdir bin
	gcc -g -Ofast -march=native -flto -Wall $(SOURCES) -lpthread -o bin/dc
	gcc -g -Ofast -march=native -flto -Wall $(SOURCES_OPT_CHECKER) -lpthread -o bin/opt_checker
