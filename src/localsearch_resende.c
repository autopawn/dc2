#include "localsearch.h"

#define NO_MOVEMENT (-2)

// ============================================================================
// Resende's and werneck local search functions

void update_structures(
        const rundata *run, const solution *sol, int u,
        const int *phi1, const int *phi2, const availmoves *avail,
        double *loss, double *gain, fastmat *extra, int undo){
    //
    const problem *prob = run->prob;

    int fr    = sol->assigns[u];
    double d_phi1 = -problem_assig_value(prob,phi1[u],u);
    double d_phi2 = -problem_assig_value(prob,phi2[u],u);
    assert(fr>=0 && avail->used[fr]);
    assert(d_phi2>=d_phi1);
    assert(phi1[u]==fr);

    if(avail->avail_rems[fr]){
        if(!undo){
            loss[fr] += d_phi2 - d_phi1;
        }else{
            loss[fr] -= d_phi2 - d_phi1;
            assert(loss[fr]>=-1e-6);
        }
    }

    // Choose between numerating possible insertions using proximity order or available insertions
    int proximity_mode = prob->n_facs/sol->n_facs <= avail->n_insertions;


    for(int k=0;k<prob->n_facs;k++){
        int fi;
        double d_fi;

        if(proximity_mode){
            fi = run->nearly_indexes[u][k];
            if(!avail->avail_inss[fi]) continue;

            d_fi = -problem_assig_value(prob,fi,u);
            if(d_fi >= d_phi2) break;
        }else{
            if(k >= avail->n_insertions) break;

            fi = avail->insertions[k];
            assert(avail->avail_inss[fi]);

            d_fi = -problem_assig_value(prob,fi,u);
            if(d_fi >= d_phi2) continue;
        }


        if(d_fi < d_phi1){
            if(!undo){
                gain[fi] += d_phi1 - d_fi;
                if(avail->avail_rems[fr]) fastmat_add(extra,fi,fr,d_phi2-d_phi1);
            }else{
                gain[fi] -= d_phi1 - d_fi;
                assert(gain[fi]>=-1e-6);
                if(avail->avail_rems[fr]) fastmat_rem(extra,fi,fr,d_phi2-d_phi1);
            }
        }else{
            if(!undo){
                if(avail->avail_rems[fr]) fastmat_add(extra,fi,fr,d_phi2-d_fi);
            }else{
                if(avail->avail_rems[fr]) fastmat_rem(extra,fi,fr,d_phi2-d_fi);
            }
        }
    }
}

double find_best_neighboor(
        const rundata *run,
        const double *loss, const double *gain, fastmat *extra, const availmoves *avail,
        int size_increase, int size_decrease, int *out_fins, int *out_frem){
    //
    const problem *prob = run->prob;
    //
    int best_fins = NO_MOVEMENT;
    int best_frem = NO_MOVEMENT;
    double best_delta = (avail->path_relinking)? -INFINITY :  0;
    // Consider just simple insertions without deletions
    if(size_increase){
        for(int t=0;t<avail->n_insertions;t++){
            int fi = avail->insertions[t];
            double delta = gain[fi] - prob->facility_cost[fi];
            if(delta > best_delta){
                best_fins  = fi;
                best_frem  = -1;
                best_delta = delta;
            }
        }
    }
    // Consider just a simple removal without insertions
    if(size_decrease){
        for(int t=0;t<avail->n_removals;t++){
            int fr = avail->removals[t];
            double delta = prob->facility_cost[fr] - loss[fr];
            if(delta > best_delta){
                best_fins = -1;
                best_frem = fr;
                best_delta = delta;
            }
        }
    }
    // Consider swaps
    coordcelltable_iter it = coordcelltable_begin(extra->tab);
    while(!coordcelltable_iter_done(&it)){
        // Get key and value from next iteration
        coord co = coordcelltable_iter_key(&it);
        cell ce = coordcelltable_iter_val(&it);
        double extrav = ce.value;
        //
        int fi = co.y;
        int fr = co.x;
        assert(ce.n_clis>0);
        assert(avail->used[fr]);
        assert(!avail->used[fi]);
        double delta = gain[fi] - loss[fr] + extrav - prob->facility_cost[fi] + prob->facility_cost[fr];
        if(delta > best_delta){
            best_fins = fi;
            best_frem = fr;
            best_delta = delta;
        }

        // Step to next iteration
        coordcelltable_iter_next(&it);
    }
    // Retrieve results
    *out_fins = best_fins;
    *out_frem = best_frem;
    return best_delta;
}

int solution_resendewerneck_hill_climbing(const rundata *run, solution **solp,const solution *target, fastmat *zeroini_mat){
    solution *sol = *solp;
    const problem *prob = run->prob;
    if(sol->n_facs<2) return 0;
    // First and Second nearest facility to each client
    int *phi1 = safe_malloc(sizeof(int)*prob->n_clis);
    int *phi2 = safe_malloc(sizeof(int)*prob->n_clis);
    for(int i=0;i<prob->n_clis;i++){
        phi1[i] = sol->assigns[i];
        phi2[i] = solution_client_2nd_nearest(prob,sol,i);
    }

    // Structures
    fastmat *extra = zeroini_mat;
    double *gain = safe_malloc(sizeof(double)*prob->n_facs);
    double *loss = safe_malloc(sizeof(double)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        gain[i] = 0;
        loss[i] = 0;
    }

    // Available moves
    availmoves *avail = availmoves_init(prob,sol,target);

    // Array of affected clients
    int n_affected = prob->n_clis;
    int *affected = safe_malloc(sizeof(int)*prob->n_clis);
    for(int i=0;i<n_affected;i++){
        affected[i] = i;
    }

    // Each movement:
    int n_moves = 0;
    int best_rem, best_ins;
    double best_delta = 0;

    int *affected_mask = safe_malloc(sizeof(int)*prob->n_clis);

    // Best solution found so far (if doing path relinking):
    solution *best_sol = NULL;
    if(avail->path_relinking){
        best_sol = solution_copy(prob,sol);
    }

    while(1){
        best_rem = NO_MOVEMENT;
        best_ins = NO_MOVEMENT;

        // Update structures for affected clients
        for(int i=0;i<n_affected;i++){
            int u = affected[i];
            update_structures(run,sol,u,phi1,phi2,avail,loss,gain,extra,0);
        }

        // Check movements that are allowed and not allowed
        int allow_size_decrease = run->local_search_rem_movement && sol->n_facs>2 &&
            (prob->size_restriction_minimum==-1 || sol->n_facs>prob->size_restriction_minimum);
        int allow_size_increase = run->local_search_add_movement &&
            (prob->size_restriction_maximum==-1 || sol->n_facs<prob->size_restriction_maximum);
        if(target){
            allow_size_increase = 1;
            allow_size_decrease = 1;
        }

        best_delta = find_best_neighboor(run,loss,gain,extra,avail,
                allow_size_increase,allow_size_decrease,&best_ins,&best_rem);

        // Stop when no movement results in a better solution:
        if(best_ins==NO_MOVEMENT){
            assert(best_rem==NO_MOVEMENT);
            break;
        }

        assert(best_rem<0 || avail->used[best_rem]);
        // Update array of affected users
        n_affected = 0;
        for(int u=0;u<prob->n_clis;u++){
            assert(sol->assigns[u] == phi1[u]);
            affected_mask[u] = 0;
            if(phi1[u]==best_rem || phi2[u]==best_rem || (-problem_assig_value(prob,best_ins,u)) < (-problem_assig_value(prob,phi2[u],u))){
                affected[n_affected] = u;
                n_affected += 1;
                affected_mask[u] = 1;
            }
        }

        for(int i=0;i<n_affected;i++){
            int u = affected[i];
            assert(affected_mask[u]);
            // Undo update structures
            update_structures(run,sol,u,phi1,phi2,avail,loss,gain,extra,1);
        }

        // Perform swap:
        assert(best_ins!=NO_MOVEMENT);
        assert(best_rem!=NO_MOVEMENT);

        double old_value = sol->value;
        if(best_rem!=-1){
            solution_remove(prob,sol,best_rem,phi2,affected_mask);
        }
        if(best_ins!=-1){
            solution_add(prob,sol,best_ins,affected_mask);
        }
        availmoves_register_move(avail,best_ins,best_rem);

        // Update best solution so far (if doing path relinking)
        if(avail->path_relinking){
            assert(best_sol);
            if(best_sol->value < sol->value){
                solution_free(best_sol);
                best_sol = solution_copy(prob,sol);
            }
        }

        assert(avail->path_relinking || sol->value>old_value);
        #ifdef DEBUG
            if(best_delta!=-INFINITY){
                double error_pred = (sol->value-old_value)-best_delta;
                if(error_pred<0) error_pred *= -1;
                assert(error_pred<1e-5 || isnan(error_pred));
            }
        #endif

        // Update phi1 and phi2
        update_phi1_and_phi2(prob,sol,best_ins,best_rem,phi1,phi2,affected_mask);

        // Count one move:
        n_moves += 1;
    }
    assert(best_delta==0 || best_delta==-INFINITY);

    free(affected_mask);

    // Clean the fastmat for reuse
    fastmat_clean(extra);

    assert((avail->path_relinking!=0) != (best_sol==NULL));

    free(affected);
    availmoves_free(avail);
    free(loss);
    free(gain);
    free(phi2);
    free(phi1);

    // Set sol to best_sol in case we are doing path relinking
    if(best_sol!=NULL){
        *solp = best_sol;
        solution_free(sol);
    }

    return n_moves;
}
