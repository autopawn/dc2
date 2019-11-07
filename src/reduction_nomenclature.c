#include "reduction.h"

const char *default_redstrategies[] = {"rand:4000","sdcebest:200"};

// Allocates an array of redstrategis from nomenclatures 
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
    int n_parts = 1;
    for(int i=0;i<nomlen;i++){
        if(nomenclature2[i]==':'){
            nomenclature2[i]=' ';
            n_parts += 1;
        }
    }

    // Compare the abreviation of the first strategy:
    char abrev[400];
    if(n_parts==2){
        if(sscanf(nomenclature2,"%s %d",abrev,&strategy.n_target)!=2){
        }
        strategy.arg = -1;
    }else if(n_parts==3){
        if(sscanf(nomenclature2,"%s %d %d",abrev,&strategy.n_target,&strategy.arg)!=3){
        }
    }else{
        fprintf(stderr,"ERROR: Invalid nomenclature format!\n");
        exit(1);
    }

    // Set the strategy method based on the abrev
    if(strcmp(abrev,"best")==0){
        strategy.method = REDUCTION_BESTS;
        assert(strategy.arg==-1);
    }
    else if(strcmp(abrev,"rand")==0){
        strategy.method = REDUCTION_RANDOM_UNIFORM;
        assert(strategy.arg==-1);
    }
    else if(strcmp(abrev,"rank")==0){
        strategy.method = REDUCTION_RANDOM_RANK;
        assert(strategy.arg==-1);
    }
    else if(strcmp(abrev,"vrh")==0){
        strategy.method = REDUCTION_VRHEURISTIC;
        if(strategy.arg==-1) strategy.arg = 2*strategy.n_target;
    }
    else if(strcmp(abrev,"sdce")==0){
        strategy.method = REDUCTION_GLOVER_SDCE;
        assert(strategy.arg==-1);
    }
    else if(strcmp(abrev,"sdcebest")==0){
        strategy.method = REDUCTION_GLOVER_SDCE_BESTS;
        assert(strategy.arg==-1);
    }
    else{
        fprintf(stderr,"ERROR: Invalid nomenclature abrev \"%s\"!\n",abrev);
        exit(1);
    }

    return strategy;
}

void redstrategy_reduce(const problem *prob, const redstrategy rstrat, 
        solution **sols, int *n_sols){
    
    printf("Reducing %5d -> %5d solutions, ",*n_sols,rstrat.n_target);
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