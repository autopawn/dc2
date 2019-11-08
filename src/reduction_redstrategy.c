#include "reduction.h"

const char *default_redstrategies[] = {"rand:4000","sdce+:200"};
 
redstrategy *redstrategy_init_from_nomenclatures(const char **noms, int *n_noms){
    int n_strategies = *n_noms;
    redstrategy *strategies; 
    if(n_strategies==0){
        // Default strategies
        n_strategies = sizeof(default_redstrategies)/sizeof(default_redstrategies[0]);
        strategies = safe_malloc(sizeof(redstrategy)*n_strategies);
        for(int i=0;i<n_strategies;i++){
            strategies[i] = redstrategy_from_nomenclature(default_redstrategies[i]);
        }
    }else{
        // Strategies from arguments
        strategies = safe_malloc(sizeof(redstrategy)*n_strategies);
        for(int i=0;i<n_strategies;i++){
            strategies[i] = redstrategy_from_nomenclature(noms[i]);
        }
    }
    // Check that strategy n_targets are decreasing
    for(int i=1;i<n_strategies;i++){
        if(strategies[i-1].n_target<=strategies[i].n_target){
            fprintf(stderr,"ERROR: strategy n_targets must decrease!\n");
        }
    }
    *n_noms = n_strategies;
    return strategies;
}

redstrategy redstrategy_from_nomenclature(const char *nomenclature){
    redstrategy strategy;
    strategy.nomenclature = nomenclature;

    int nomlen = strlen(nomenclature);
    assert(nomlen<400);
    char nomenclature2[400];
    strcpy(nomenclature2,nomenclature);

    // Replace the colons with spaces
    int strategy_n_parts = 1;
    for(int i=0;i<nomlen;i++){
        if(nomenclature2[i]==':'){
            nomenclature2[i]=' ';
            strategy_n_parts += 1;
        }
    }

    // Get the abreviations of the strategies:
    char abrev[400];
    char distm[400];
    
    int n_scan = sscanf(nomenclature2,"%s %d %s %d",
            abrev,&strategy.n_target,distm,&strategy.arg);
    if(n_scan<2){
        fprintf(stderr,"ERROR: Invalid reduction format \"%s\"!\n",nomenclature);
        exit(1);
    }

    // Set the strategy arg to -1 if not given:
    if(n_scan<4) strategy.arg = -1;

    int invalid = 0;
    // Set the strategy method based on the abrev
    if(strcmp(abrev,"best")==0){
        strategy.method = REDUCTION_BESTS;
        if(strategy_n_parts!=2) invalid = 1;
    }
    else if(strcmp(abrev,"rand")==0){
        strategy.method = REDUCTION_RANDOM_UNIFORM;
        if(strategy_n_parts!=2) invalid = 1;
    }
    else if(strcmp(abrev,"rank")==0){
        strategy.method = REDUCTION_RANDOM_RANK;
        if(strategy_n_parts!=2) invalid = 1;
    }
    else if(strcmp(abrev,"vrh")==0){
        strategy.method = REDUCTION_VRHEURISTIC;
        if(strategy_n_parts>4) invalid = 1;
        // Set the default value for the VISION RANGE
        if(strategy.arg==-1) strategy.arg = 2*strategy.n_target;
    }
    else if(strcmp(abrev,"sdce")==0){
        strategy.method = REDUCTION_GLOVER_SDCE;
        if(strategy_n_parts!=3) invalid = 1;
    }
    else if(strcmp(abrev,"sdce+")==0){
        strategy.method = REDUCTION_GLOVER_SDCE_BESTS;
        if(strategy_n_parts!=3) invalid = 1;
    }
    else{
        fprintf(stderr,"ERROR: Invalid reduction name \"%s\"!\n",abrev);
        exit(1);
    }
    if(invalid){
        fprintf(stderr,"ERROR: Invalid reduction arguments \"%s\"!\n",nomenclature);
        exit(1);
    }

    // Identify the dissimilitude and distance strategies
    if(n_scan<3){
        strategy.soldis = SOLDIS_MEAN_SQUARE_ERROR;
        strategy.facdis = FACDIS_SUM_OF_DELTAS;
    }else{
        if(strcmp(distm,"msemin")==0){
            strategy.soldis = SOLDIS_MEAN_SQUARE_ERROR;
            strategy.facdis = FACDIS_MIN_TRIANGLE;
        }
        else if(strcmp(distm,"msesum")==0){
            strategy.soldis = SOLDIS_MEAN_SQUARE_ERROR;
            strategy.facdis = FACDIS_SUM_OF_DELTAS;
        }
        else if(strcmp(distm,"haumin")==0){
            strategy.soldis = SOLDIS_HAUSDORF;
            strategy.facdis = FACDIS_MIN_TRIANGLE;
        }
        else if(strcmp(distm,"hausum")==0){
            strategy.soldis = SOLDIS_HAUSDORF;
            strategy.facdis = FACDIS_SUM_OF_DELTAS;
        }
        else if(strcmp(distm,"pcd")==0){
            strategy.soldis = SOLDIS_PER_CLIENT_DELTA;
            strategy.facdis = 0; // Placeholder.
        }
        else{
            fprintf(stderr,"ERROR: Invalid dissimilitude abrev \"%s\"!\n",distm);
            exit(1);
        }
    }

    return strategy;
}

void redstrategy_reduce(const problem *prob, const redstrategy rstrat, 
        solution **sols, int *n_sols){

    if(*n_sols<=rstrat.n_target) return;
    
    printf("Reducing \033[31;1m%d\033[0m -> \033[31;1m%d\033[0m solutions, ",*n_sols,rstrat.n_target);
    if(rstrat.method==REDUCTION_BESTS){
        printf("selecting bests.\n");
        reduction_bests(prob,sols,n_sols,rstrat.n_target);
    }
    else if(rstrat.method==REDUCTION_RANDOM_UNIFORM){
        printf("randomly (uniform).\n");
        reduction_random_uniform(prob,sols,n_sols,rstrat.n_target);
    }
    else if(rstrat.method==REDUCTION_RANDOM_RANK){
        printf("randomly (by rank).\n");
        reduction_random_rank(prob,sols,n_sols,rstrat.n_target);
    }
    else{
        printf("???.\n");
        fprintf(stderr,"ERROR: Reduction method not yet implemented.\n");
        exit(1);
    }
    
}