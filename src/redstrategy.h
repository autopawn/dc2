#ifndef DC_REDSTRATEGY_H
#define DC_REDSTRATEGY_H

#include "utils.h"

#define N_FACDIS_MODES 3

typedef enum {
    FACDIS_MIN_TRIANGLE = 0,  // df(a,b) = min_j d(a,j)+d(b,j)
    FACDIS_SUM_OF_DELTAS = 1, // df(a,b) = sum_j |d(a,j)-d(b,j)|
    FACDIS_NONE = 2,
} facdismode;

typedef enum {
    SOLDIS_MEAN_SQUARE_ERROR = 0, // D(A,B) = sum_j |d(A,j)-d(B,j)|
    SOLDIS_HAUSDORF = 1,          // D(A,B) = max {sup_a inf_b df(a,b), sup_b inf_a df(b,a)}
    SOLDIS_PER_CLIENT_DELTA = 2,  // D(A,B) = sum_a inf_b df(a,b) + sum_b inf_a df(b,a)
} soldismode;

typedef enum {
    // | Pick the best solutions
    REDUCTION_BESTS,
    // | Random uniform
    REDUCTION_RANDOM_UNIFORM,
    // | Random with ponderation by rank
    REDUCTION_RANDOM_RANK,
    // | The VR-Heuristic
    REDUCTION_VRHEURISTIC,
    // | Enhanced Glover simple diversity-based starting method
    REDUCTION_GLOVER_SDCE,
    // | EGLOVER_SDCE selecting the best solution of each cluster 
    REDUCTION_GLOVER_SDCE_BESTS,
} reduction_method;

typedef struct {
    const char *nomenclature;
    reduction_method method;
    int n_target;
    soldismode soldis;
    facdismode facdis;
    int arg;
} redstrategy;

// Allocates an array of redstrategies from nomenclatures 
redstrategy *redstrategy_init_from_nomenclatures(const char **noms, int *n_noms);

// Parses a nomenclature to generate a redstrategy
redstrategy redstrategy_from_nomenclature(const char *nomenclature);

// Gets the facdis mode required by this redstrategy
facdismode redstrategy_required_facdis_mode(const redstrategy rstrat);

#endif