#include "solution.h"

// Used to get the size of the linear-growing array reserved for facility indexes.
int get_size(int len){
    const int FACTOR = 8;
    int s = FACTOR*((len+FACTOR-1)/FACTOR);
    return s==0? FACTOR : s;
}

// Compare solutions to sort on decreasing value
int solutionp_value_cmp_inv(const void *a, const void *b){
    solution **aa = (solution **)a;
    solution **bb = (solution **)b;
    double diff = (*bb)->value - (*aa)->value;
    return diff == 0 ? 0 : ((diff > 0) ? +1 : -1);
}

solution *solution_empty(const problem *prob){
    solution *sol = safe_malloc(sizeof(solution));
    sol->n_facs = 0;
    sol->facs = safe_malloc(sizeof(int)*get_size(prob->n_facs));
    sol->assigns = safe_malloc(sizeof(int)*prob->n_clis);
    for(int j=0;j<prob->n_clis;j++){
        sol->assigns[j] = -1;
    }
    // Initialize solution value
    sol->value = 0;
    for(int j=0;j<prob->n_clis;j++){
        sol->value += problem_assig_value(prob,-1,j);
    }
    return sol;
}

solution *solution_copy(const problem *prob, const solution *sol){
    solution *sol2 = safe_malloc(sizeof(solution));
    sol2->n_facs = sol->n_facs;
    sol2->facs = safe_malloc(sizeof(int)*get_size(prob->n_facs));
    memcpy(sol2->facs,sol->facs,sizeof(int)*sol->n_facs);
    sol2->assigns = safe_malloc(sizeof(int)*prob->n_clis);
    memcpy(sol2->assigns,sol->assigns,sizeof(int)*prob->n_clis);
    sol2->value = sol->value;
    return sol2;
}

void solution_add(const problem *prob, solution *sol, int newf){
    // Check if f is already on the solution:
    for(int f=0;f<sol->n_facs;f++){
        if(sol->facs[f]==newf) return;
    }
    // Extend array of facility indexes if necessary
    if(sol->n_facs==get_size(sol->n_facs)){
        sol->facs = realloc(sol->facs,sizeof(int)*get_size(sol->n_facs+1));
    }
    // Add facility to the solution
    add_to_sorted(sol->facs,&sol->n_facs,newf);
    // | New value after adding the new facility.
    double value2 = 0;
    // Reassign clients to the new instalation
    for(int c=0;c<prob->n_clis;c++){
        double val_pre = problem_assig_value(prob,sol->assigns[c],c);
        double val_pos = problem_assig_value(prob,newf,c);
        if(val_pos>val_pre){
            sol->assigns[c] = newf;
            value2 += val_pos;
        }else{
            value2 += val_pre;
        }
    }
    // Cost of facilities
    for(int i=0;i<sol->n_facs;i++){
        value2 -= prob->facility_cost[i];
    }
    // Update solution value
    sol->value = value2;
}

void solution_free(solution *sol){
    free(sol->facs);
    free(sol->assigns);
    free(sol);
}