#include "solution.h"

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
    sol->facs = safe_malloc(sizeof(int)*1);
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
    sol2->facs =    safe_malloc(sizeof(int)*sol->n_facs);
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
    // Extend array of facilities
    sol->facs = realloc(sol->facs,sizeof(int)*(sol->n_facs+1));
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

void solution_remove(const problem *prob, solution *sol, int remf){
    rem_of_sorted(sol->facs,&sol->n_facs,remf);
    // New value after adding the new facility
    double value2 = 0;
    // Drop clients of the facility.
    for(int c=0;c<prob->n_clis;c++){
        // If the client was owned by the facility reassing
        if(sol->assigns[c]==remf){
            int reassign = -1;
            double reassign_value = problem_assig_value(prob,reassign,c);
            for(int i=0;i<sol->n_facs;i++){
                int candidate = sol->facs[i];
                double cand_value = problem_assig_value(prob,candidate,c);
                if(cand_value>reassign_value){
                    reassign = candidate;
                    reassign_value = cand_value;
                }
            }
            // Reassign client
            sol->assigns[c] = reassign;
            // Value of new assigment
            value2 += reassign_value;
        }else{
            // Value of the current assignment
            value2 += problem_assig_value(prob,sol->assigns[c],c);
        }
    }
    // The facility costs
    for(int i=0;i<sol->n_facs;i++){
        value2 -= prob->facility_cost[sol->facs[i]];
    }
    // Update solution value
    sol->value = value2;
}

void solution_free(solution *sol){
    free(sol->facs);
    free(sol->assigns);
    free(sol);
}

// Compute the distance between two solutions
double solution_dissimilitude(const problem *prob,
        const solution *sol1, const solution *sol2,
        soldismode sdismode, facdismode fdismode){
    // Compute the dissimilitude according to the sdismode
    if(sdismode==SOLDIS_MEAN_SQUARE_ERROR){
        // Expect the facility distances for this mode to be computed:
        if(prob->facs_distance[fdismode]==NULL){
            fprintf(stderr,"Error: problem facility-facility distances are not precomputed!\n");
            exit(1);
        }
        // Add distance from each facility in sol1 to sol2.
        double disim = 0;
        for(int t=0;t<2;t++){
            for(int i1=0;i1<sol1->n_facs;i1++){
                double min_dist = INFINITY;
                int f1 = sol1->facs[i1];
                for(int i2=0;i2<sol2->n_facs;i2++){
                    int f2 = sol2->facs[i2];
                    double dist = prob->facs_distance[fdismode][f1][f2];
                    if(dist<min_dist) min_dist = dist;
                }
                disim += min_dist;
            }
            // Swap solutions for 2nd iteration:
            const solution *aux = sol1; sol1 = sol2; sol2 = aux;
        }
        return disim;
    }
    else if(sdismode==SOLDIS_HAUSDORF){ // Based on https://github.com/mavillan/py-hausdorff
        // Expect the facility distances for this mode to be computed:
        if(prob->facs_distance[fdismode]==NULL){
            fprintf(stderr,"Error: problem facility-facility distances are not precomputed!\n");
            exit(1);
        }
        double disim = 0;
        for(int t=0;t<2;t++){
            for(int i1=0;i1<sol1->n_facs;i1++){
                int f1 = sol1->facs[i1];
                double cmin = INFINITY;
                for(int i2=0;i2<sol2->n_facs;i2++){
                    int f2 = sol2->facs[i2];
                    double dist = prob->facs_distance[fdismode][f1][f2];
                    if(dist<cmin) cmin = dist;
                    if(cmin<disim) break;
                }
                if(disim<cmin && cmin<INFINITY) disim = cmin;
            }
            // Swap solutions for 2nd iteration:
            const solution *aux = sol1; sol1 = sol2; sol2 = aux;
        }
        return disim;
    }
    else if(sdismode==SOLDIS_PER_CLIENT_DELTA){
        double total = 0;
        for(int i=0;i<prob->n_clis;i++){
            double cost_a = prob->distance[sol1->assigns[i]][i];
            double cost_b = prob->distance[sol2->assigns[i]][i];
            double delta = cost_a-cost_b;
            if(delta<0) delta = -delta;
            total += delta;
        }
        return total;
    }
    // Error if it was an invalid soldismode
    fprintf(stderr,"ERROR: Invalid soldismode!\n");
    exit(1);
}

// An upper bound for the best value that a children solution could have
double solution_upper_bound(const problem *prob, const solution *sol){
    double upbound = prob->precomp_client_optimal_gain;
    for(int i=0;i<sol->n_facs;i++){
        upbound -= prob->facility_cost[sol->facs[i]];
    }
    return upbound;
}

void solution_print(const problem *prob, const solution *sol, FILE *fp){
    fprintf(fp,"== SOLUTION ==\n");
    fprintf(fp,"# VALUE: %lf\n",sol->value);
    fprintf(fp,"# ASSIGNS:");
    for(int i=0;i<prob->n_clis;i++){
        fprintf(fp," %d",sol->assigns[i]);
    }
    fprintf(fp,"\n");
    fprintf(fp,"# N_FACS: %d\n",sol->n_facs);
    //
    fprintf(fp,"# INDEXES: ");
    for(int i=0;i<sol->n_facs;i++){
        fprintf(fp,"%d ",sol->facs[i]);
    }
    fprintf(fp,"\n");
    // Print clients for each facility
    for(int i=0;i<sol->n_facs;i++){
        fprintf(fp,"FAC %d :",sol->facs[i]);
        for(int j=0;j<prob->n_clis;j++){
            if(sol->assigns[j]==sol->facs[i]){
                fprintf(fp," %d",j);
            }
        }
        fprintf(fp,"\n");
    }
}