#include "redstrategy.h"

const char *default_redstrategies[] = {"rand1:3000","sdbs+:100:pcd","_best:200"};

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
    for(int t=0;t<=1;t++){ // For not PR and for PR
        int current_n_target = INT_MAX;
        for(int i=0;i<n_strategies;i++){
            if(strategies[i].for_selected_sols == t){
                if(strategies[i].n_target>=current_n_target){
                    fprintf(stderr,"ERROR: strategy n_targets must decrease!\n");
                    exit(1);
                }
                current_n_target = strategies[i].n_target;
            }
        }
    }

    *n_noms = n_strategies;
    return strategies;
}

redstrategy redstrategy_from_nomenclature(const char *nomenclature){
    redstrategy strategy;
    strategy.nomenclature = nomenclature;

    // Check if reduction method starts with '_', meaning that it is intended for path relinking
    if(nomenclature[0]=='_'){
        strategy.for_selected_sols = 1;
        nomenclature = &nomenclature[1];
    }else{
        strategy.for_selected_sols = 0;
    }

    //

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
        fprintf(stderr,"ERROR: Invalid reduction format \"%s\"!\n",strategy.nomenclature);
        exit(1);
    }

    // Set the strategy arg to -1 if not given:
    if(n_scan<4) strategy.arg = -1;
    // Not elitist by default
    strategy.elitist = 0;

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
    else if(strcmp(abrev,"rand1")==0){
        strategy.method = REDUCTION_RANDOM_UNIFORM;
        if(strategy_n_parts!=2) invalid = 1;
        strategy.elitist = 1;
    }
    else if(strcmp(abrev,"rank")==0){
        strategy.method = REDUCTION_RANDOM_RANK;
        if(strategy_n_parts!=2) invalid = 1;
    }
    else if(strcmp(abrev,"rank1")==0){
        strategy.method = REDUCTION_RANDOM_RANK;
        if(strategy_n_parts!=2) invalid = 1;
        strategy.elitist = 1;
    }
    else if(strcmp(abrev,"vrh")==0){
        strategy.method = REDUCTION_VRHEURISTIC;
        if(strategy_n_parts>4) invalid = 1;
        // Set the default value for the VISION RANGE
        if(strategy.arg==-1) strategy.arg = 2*strategy.n_target;
    }
    else if(strcmp(abrev,"sdbs")==0){
        strategy.method = REDUCTION_GLOVER_SDBS;
        if(strategy_n_parts>3) invalid = 1;
    }
    else if(strcmp(abrev,"sdbs+")==0){
        strategy.method = REDUCTION_GLOVER_SDBS_BESTS;
        if(strategy_n_parts>3) invalid = 1;
    }
    else{
        fprintf(stderr,"ERROR: Invalid reduction name \"%s\"!\n",abrev);
        exit(1);
    }
    if(invalid){
        fprintf(stderr,"ERROR: Invalid reduction arguments \"%s\"!\n",strategy.nomenclature);
        exit(1);
    }

    // Identify the dissimilitude and distance strategies
    if(n_scan<3){
        // Default dissimilitude:
        strategy.soldis = SOLDIS_PER_CLIENT_DELTA;
        strategy.facdis = FACDIS_NONE;
    }else{
        if(strcmp(distm,"mgemin")==0){
            strategy.soldis = SOLDIS_MEAN_GEOMETRIC_ERROR;
            strategy.facdis = FACDIS_MIN_TRIANGLE;
        }
        else if(strcmp(distm,"mgesum")==0){
            strategy.soldis = SOLDIS_MEAN_GEOMETRIC_ERROR;
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
            strategy.facdis = FACDIS_NONE;
        }
        else if(strcmp(distm,"autosum")==0){
            strategy.soldis = SOLDIS_AUTO;
            strategy.facdis = FACDIS_SUM_OF_DELTAS;
        }
        else if(strcmp(distm,"automin")==0){
            strategy.soldis = SOLDIS_AUTO;
            strategy.facdis = FACDIS_MIN_TRIANGLE;
        }
        else if(strcmp(distm,"indexval")==0){
            strategy.soldis = SOLDIS_INDEXES_VALUE;
            strategy.facdis = FACDIS_NONE;
        }
        else{
            fprintf(stderr,"ERROR: Invalid dissimilitude abrev \"%s\"!\n",distm);
            exit(1);
        }
    }

    return strategy;
}

facdismode redstrategy_required_facdis_mode(const redstrategy rstrat){
    /* Check if the reduction method doesn't use dissimilitude */
    if(rstrat.method==REDUCTION_BESTS) return FACDIS_NONE;
    if(rstrat.method==REDUCTION_RANDOM_UNIFORM) return FACDIS_NONE;
    if(rstrat.method==REDUCTION_RANDOM_RANK) return FACDIS_NONE;
    return rstrat.facdis;
}
