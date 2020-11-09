#include "solution.h"

// Compare solutions to sort on decreasing value
int solutionp_value_cmp_inv(const void *a, const void *b){
    const solution **aa = (const solution **)a;
    const solution **bb = (const solution **)b;
    double diff = (*bb)->value - (*aa)->value;
    return diff == 0 ? 0 : ((diff > 0) ? +1 : -1);
}

// Compare solutions for equality
int solutionp_facs_cmp(const void *a, const void *b){
    const solution **aa = (const solution **) a;
    const solution **bb = (const solution **) b;
    const solution *sol1 = *aa;
    const solution *sol2 = *bb;
    int d = sol1->n_facs - sol2->n_facs;
    if(d!=0) return d;
    for(int i=0;i<sol1->n_facs;i++){
        d = sol1->facs[i]-sol2->facs[i];
        if(d!=0) return d;
    }
    return 0;
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
    sol->terminal = 0;
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
    sol2->terminal = sol->terminal;
    return sol2;
}

void solution_add(const problem *prob, solution *sol, int newf, int *affected){
    // Check if f is already on the solution:
    for(int f=0;f<sol->n_facs;f++){
        if(sol->facs[f]==newf){
            return;
        }
    }
    // Extend array of facilities
    sol->facs = safe_realloc(sol->facs,sizeof(int)*(sol->n_facs+1));
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
        value2 -= prob->facility_cost[sol->facs[i]];
    }
    // Update solution value
    sol->value = value2;
}

void solution_remove(const problem *prob, solution *sol, int remf, int *phi2, int *affected){
    rem_of_sorted(sol->facs,&sol->n_facs,remf);
    // New value after adding the new facility
    double value2 = 0;
    // Drop clients of the facility.
    for(int c=0;c<prob->n_clis;c++){
        // If the client was owned by the facility reassing
        if(sol->assigns[c]==remf){
            int reassign;
            double reassign_value;
            if(phi2!=NULL){
                reassign = phi2[c];
                reassign_value = problem_assig_value(prob,reassign,c);
            }else{
                reassign = -1;
                reassign_value = problem_assig_value(prob,reassign,c);
                for(int i=0;i<sol->n_facs;i++){
                    int candidate = sol->facs[i];
                    double cand_value = problem_assig_value(prob,candidate,c);
                    if(cand_value>reassign_value){
                        reassign = candidate;
                        reassign_value = cand_value;
                    }
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
double solution_dissimilitude(const rundata *run,
        const solution *sol1, const solution *sol2,
        soldismode sdismode, facdismode fdismode){
    // The auto case picks the sdismode
    if(sdismode==SOLDIS_AUTO){
        if(sol1->n_facs*sol2->n_facs <= 15*run->prob->n_clis){
            sdismode = SOLDIS_MEAN_GEOMETRIC_ERROR;
        }else{
            sdismode = SOLDIS_PER_CLIENT_DELTA;
        }
    }
    // Compute the dissimilitude according to the sdismode
    if(sdismode==SOLDIS_MEAN_GEOMETRIC_ERROR){
        // Expect the facility distances for this mode to be computed:
        if(run->precomp->facs_distance[fdismode]==NULL){
            fprintf(stderr,"Error: problem facility-facility distances are not precomputed!\n");
            exit(1);
        }
        // Add distance from each facility in sol1 to sol2.
        double disim = 0;
        int i1 = 0;
        int i2 = 0;
        // Make use of the fact that facilities are sorted in both solutions
        while(i1<sol1->n_facs || i2<sol2->n_facs){
            int s1f = i1<sol1->n_facs? sol1->facs[i1] : INT_MAX;
            int s2f = i2<sol2->n_facs? sol2->facs[i2] : INT_MAX;
            if(s1f==s2f){
                // No delta added, same facility
                i1 += 1;
                i2 += 1;
            }else if(s1f<s2f){
                // Add delta from s1f to sol2
                double min_dist = INFINITY;
                for(int k=0;k<sol2->n_facs;k++){
                    int f2 = sol2->facs[k];
                    double dist = run->precomp->facs_distance[fdismode][s1f][f2];
                    if(dist<min_dist) min_dist = dist;
                }
                disim += min_dist;
                //
                i1 += 1;
            }else{
                // Add delta from s2f to sol1
                double min_dist = INFINITY;
                for(int k=0;k<sol1->n_facs;k++){
                    int f1 = sol1->facs[k];
                    double dist = run->precomp->facs_distance[fdismode][s2f][f1];
                    if(dist<min_dist) min_dist = dist;
                }
                disim += min_dist;
                //
                i2 += 1;
            }
        }
        return disim/(sol1->n_facs+sol2->n_facs);
        // // OLDER VERSION
        // double disim = 0;
        // for(int t=0;t<2;t++){
        //     for(int i1=0;i1<sol1->n_facs;i1++){
        //         double min_dist = INFINITY;
        //         int f1 = sol1->facs[i1];
        //         for(int i2=0;i2<sol2->n_facs;i2++){
        //             int f2 = sol2->facs[i2];
        //             double dist = run->precomp->facs_distance[fdismode][f1][f2];
        //             if(dist<min_dist) min_dist = dist;
        //         }
        //         disim += min_dist;
        //     }
        //     // Swap solutions for 2nd iteration:
        //     const solution *aux = sol1; sol1 = sol2; sol2 = aux;
        // }
    }
    else if(sdismode==SOLDIS_HAUSDORF){ // Based on https://github.com/mavillan/py-hausdorff
        // Expect the facility distances for this mode to be computed:
        if(run->precomp->facs_distance[fdismode]==NULL){
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
                    double dist = run->precomp->facs_distance[fdismode][f1][f2];
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
        for(int i=0;i<run->prob->n_clis;i++){
            double cost_a = problem_assig_value(run->prob,sol1->assigns[i],i);
            double cost_b = problem_assig_value(run->prob,sol2->assigns[i],i);
            double delta = cost_a-cost_b;
            if(delta<0) delta = -delta;
            total += delta;
        }
        return total;
    }else if(sdismode==SOLDIS_INDEXES_VALUE){
        int delta = diff_sorted(sol1->facs,sol1->n_facs,sol2->facs,sol2->n_facs);
        double value_delta = sol1->value - sol2->value;
        if(value_delta<0) value_delta *= -1;
        // Use the number of different indexes as dissimilitude and the value to break ties
        double delta_val = 67108864*delta + value_delta;
        return delta_val;
    }
    // Error if it was an invalid soldismode
    fprintf(stderr,"ERROR: Invalid soldismode!\n");
    exit(1);
}

// An upper bound for the best value that a children solution could have
double solution_upper_bound(const rundata *run, const solution *sol){
    double upbound = run->precomp->precomp_client_optimal_gain;
    for(int i=0;i<sol->n_facs;i++){
        upbound -= run->prob->facility_cost[sol->facs[i]];
    }
    return upbound;
}

int solution_client_2nd_nearest(const problem *prob, const solution *sol, int cli){
    // Nearest facility
    int phi1 = sol->assigns[cli];
    // 2nd nearest
    int phi2 = -1;
    double phi2_assig_cost = problem_assig_cost(prob,phi2,cli);
    // Iterate to find 2nd best facility
    for(int p=0;p<sol->n_facs;p++){
        int fac = sol->facs[p];
        if(fac==phi1) continue;
        double fac_assig_cost = problem_assig_cost(prob,fac,cli);
        if(fac_assig_cost<phi2_assig_cost){
            phi2 = fac;
            phi2_assig_cost = fac_assig_cost;
        }
    }
    return phi2;
}

void solution_findout(const problem *prob, const solution *sol, int f_ins, double *v,
        const int *phi2, int *frem_allowed,
        int *out_f_rem, double *out_profit, double *out_profit_worem){
    // The gain for swaping decreases on the cost of the inserted facility
    double w = f_ins==-1? 0 : -prob->facility_cost[f_ins];
    // The cost for swapping decreases on the cost of the removed facility
    for(int k=0;k<sol->n_facs;k++){
        v[sol->facs[k]] = -prob->facility_cost[sol->facs[k]];
    }
    //
    for(int u=0;u<prob->n_clis;u++){
        int phi1u = sol->assigns[u];
        double assig_f_ins_value = problem_assig_value(prob,f_ins,u);
        double assig_phi1u_value = problem_assig_value(prob,phi1u,u);
        double delta = assig_f_ins_value - assig_phi1u_value;
        if(delta>=0){ // Profit by adding f_ins, because it is nearly.
            w += delta;
            assert(f_ins!=-1 || delta<=0.000001);
        }else{ // Loss by removing phi1u, because it is nearly.
            assert(phi1u==-1 || v[phi1u]!=-INFINITY);
            if(phi1u==-1) continue; // phi1u not part of the solution.
            double assig_phi2u_value = problem_assig_value(prob,phi2[u],u);
            if(assig_f_ins_value > assig_phi2u_value){
                v[phi1u] += assig_phi1u_value - assig_f_ins_value;
            }else{
                v[phi1u] += assig_phi1u_value - assig_phi2u_value;
            }
        }
    }
    // Find the one to be removed with less loss
    int f_rem = -1;
    for(int k=0;k<sol->n_facs;k++){
        int f = sol->facs[k];
        if(frem_allowed==NULL || frem_allowed[f]){
            if(f_rem==-1 || v[f]<v[f_rem]) f_rem = f;
        }
    }
    // Outputs
    *out_f_rem = f_rem;
    *out_profit = w - (f_rem==-1? 0 : v[f_rem]);
    *out_profit_worem = w;
    // Reset v
    for(int k=0;k<sol->n_facs;k++){
        v[sol->facs[k]] = -INFINITY;
    }
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

int solution_check_integrity(const problem *prob, const solution *sol){
    int integrity = 1;
    // Check that each client is assigned to it's nearest facility in the solution
    for(int j=0;j<prob->n_clis;j++){
        double current_value = problem_assig_value(prob,sol->assigns[j],j);
        for(int k=0;k<sol->n_facs;k++){
            int f = sol->facs[k];
            double other_value = problem_assig_value(prob,f,j);
            if(current_value < other_value) integrity = 0;
        }
    }
    // Check that he value corresponds with the stored value
    double value = 0;
    for(int j=0;j<prob->n_clis;j++){
        value += problem_assig_value(prob,sol->assigns[j],j);
    }
    for(int k=0;k<sol->n_facs;k++){
        int f = sol->facs[k];
        value -= prob->facility_cost[f];
    }
    if(isfinite(value)){
        double error = sol->value-value;
        if(error<0) error *= -1;
        if(error>=1e-5 && !isnan(error)) integrity = 0;
    }

    return integrity;
}
