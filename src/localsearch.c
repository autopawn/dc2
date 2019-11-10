#include "localsearch.h"

int solution_client_2nd_nearest(const problem *prob, const solution *sol, int cli){
    // Nearest facility
    int phi1 = sol->assigns[cli];
    // 2nd nearest
    int phi2 = -1;
    double phi2_assig_value = problem_assig_value(prob,phi2,cli);
    // Iterate to find 2nd best facility
    for(int p=0;p<sol->n_facs;p++){
        int fac = sol->facs[p];
        if(fac==phi1) continue;
        if(problem_assig_value(prob,fac,cli)>phi2_assig_value){
            phi2 = fac;
            phi2_assig_value = problem_assig_value(prob,phi2,cli);
        }
    }
    return phi2;
}

// Find the best option for removal if f_ins is inserted to the solution
// NOTE: v must be intialized with -INFINITY and have size equal to prob->n_facs.
// ^ It is always reset to that state before returning.
void solution_findout(const problem *prob, const solution *sol, int f_ins, double *v,
        const int *phi2, int *out_f_rem, double *out_profit){
    // The facility costs:
    double w = -prob->facility_cost[f_ins];
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
        }else{ // Loss by removing phi1u, because it is nearly.
            if(phi1u==-1 || v[phi1u]==-INFINITY) continue; // phi1u not part of the solution.
            double assig_phi2u_value = problem_assig_value(prob,phi2[u],u);
            if(assig_f_ins_value > assig_phi2u_value){
                v[phi1u] += assig_phi1u_value - assig_f_ins_value;
            }else{
                v[phi1u] += assig_phi1u_value - assig_phi2u_value;
            }
        }
    }
    // Find the one to be removed with less loss
    int f_rem = sol->facs[0];
    for(int k=1;k<sol->n_facs;k++){
        if(v[sol->facs[k]]<v[f_rem]){
            f_rem = sol->facs[k];
        }
    }
    // Outputs
    *out_f_rem = f_rem;
    *out_profit = w - (f_rem==-1? 0 : v[f_rem]);
    // Reset v
    for(int k=0;k<sol->n_facs;k++){
        v[sol->facs[k]] = -INFINITY;
    }
}

void solution_whitaker_hill_climbing(const problem *prob, solution *sol){
    if(sol->n_facs==0) return;
    // Second nearest facility to each client
    int *phi2 = safe_malloc(sizeof(int)*prob->n_clis);
    // Array if a facility appears in the solution
    int *used = safe_malloc(sizeof(int)*prob->n_facs);
    // Auxiliar array for solution_findout
    double *v = safe_malloc(sizeof(double)*prob->n_facs);
    // Initialize arrays
    for(int i=0;i<prob->n_facs;i++){
        used[i] = 0;
        v[i] = -INFINITY;
    }
    // Each movement:
    while(1){
        for(int i=0;i<sol->n_facs;i++) used[sol->facs[i]] = 1;
        // Clients nearest to the solution
        for(int i=0;i<prob->n_clis;i++){
            phi2[i] = solution_client_2nd_nearest(prob,sol,i);
        }
        // Insertion candidate:
        double best_delta = 0;
        int best_rem = -1;
        int best_ins = -1;
        for(int f=0;f<prob->n_facs;f++){
            if(used[f]) continue;
            // Find the best option for removal:
            int f_rem;
            double delta_profit;
            solution_findout(prob,sol,f,v,phi2,&f_rem,&delta_profit);
            if(delta_profit>best_delta){
                best_delta = delta_profit;
                best_rem = f_rem;
                best_ins = f;
            }
        }
        // Stop when no movement results in a better solution:
        if(best_ins==-1) break;
        // Perform swap:
        assert(best_rem!=-1);
        double old_value = sol->value;
        solution_remove(prob,sol,best_rem);
        used[best_rem] = 0;
        solution_add(prob,sol,best_ins);
        used[best_rem] = 1;
        assert(sol->value>old_value);
    }


    // Free memory
    free(v);
    free(used);
    free(phi2);
}
