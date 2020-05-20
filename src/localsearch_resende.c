#include "localsearch.h"

#define NO_MOVEMENT (-2)

// ============================================================================
// fastmatrix

/* A fastmatrix is a matrix that allows a fast iteration over non-zero elements. */

// coordinate in the fastmat
typedef struct {
    int y,x;
} coord;

// cell in the fastmat
typedef struct {
    // current cell value
    double value;
    // number of clients adding to this value (when 0, the value should be 0)
    int n_clis;
} cell;

// the fastmat
struct fastmat {
    // Values of the matrix's cells
    int size_y, size_x;
    cell **cells;
    // Array of non-zero entries
    int n_nonzeros;
    int s_nonzeros;
    coord *nonzeros;
};

// Initializes fastmat
fastmat *fastmat_init(int size_y, int size_x){
    fastmat *mat = safe_malloc(sizeof(fastmat));
    mat->cells = safe_malloc(sizeof(cell*)*size_y);
    for(int y=0;y<size_y;y++){
        mat->cells[y] = safe_malloc(sizeof(cell)*size_x);
        for(int x=0;x<size_x;x++){
            mat->cells[y][x].value  = 0;
            mat->cells[y][x].n_clis = 0;
        }
    }
    mat->size_y = size_y;
    mat->size_x = size_x;
    mat->s_nonzeros = size_y*size_x;
    mat->nonzeros = safe_malloc(sizeof(coord)*mat->s_nonzeros);
    mat->n_nonzeros = 0;
    return mat;
}

// Delete coords that point to a cell with value 0 from the nonzeros array.
// Should be called after removing values.
void fastmat_pack(fastmat *mat){ // Note: It doesn't merge repeated apparances (but they shall not appear).
    int new_n_nonzeros = 0;
    for(int i=0;i<mat->n_nonzeros;i++){
        coord co = mat->nonzeros[i];
        assert(mat->cells[co.y][co.x].n_clis>=0);
        if(mat->cells[co.y][co.x].n_clis>0){
            mat->nonzeros[new_n_nonzeros] = mat->nonzeros[i];
            new_n_nonzeros += 1;
        }
    }
    mat->n_nonzeros = new_n_nonzeros;
}

// Adds a value on the given position in the fastmatrix
void fastmat_add(fastmat *mat, int y, int x, double v){
    // Find current cell
    assert(y>=0 && x>=0 && y<mat->size_y && x<mat->size_x);
    assert(v>=1e-6);
    mat->cells[y][x].value += v;
    mat->cells[y][x].n_clis += 1;
    if(mat->cells[y][x].n_clis==1){ // Add nonzero to the list
        assert(mat->n_nonzeros<mat->s_nonzeros);
        coord *co = &mat->nonzeros[mat->n_nonzeros];
        co->y = y;
        co->x = x;
        mat->n_nonzeros++;
    }
}

// Intended for reverting an addition on the given position in a fastmatrix
void fastmat_rem(fastmat *mat, int y, int x, double v){
    // if(mat->cells[y][x].n_clis==0) return; // ???
    assert(mat->cells[y][x].n_clis>0);
    mat->cells[y][x].value  -= v;
    mat->cells[y][x].n_clis -= 1;
    if(mat->cells[y][x].n_clis==0){
        #ifdef DEBUG
            assert(mat->cells[y][x].value<=0.000001);
        #endif
        // This may be required due rounding errors
        mat->cells[y][x].value = 0;
    }
}

// Free fastmat memory
void fastmat_free(fastmat *mat){
    for(int y=0;y<mat->size_y;y++) free(mat->cells[y]);
    free(mat->cells);
    free(mat->nonzeros);
    free(mat);
}

void fastmat_clean(fastmat *mat){
    for(int i=0;i<mat->n_nonzeros;i++){
        coord co = mat->nonzeros[i];
        mat->cells[co.y][co.x].value  = 0;
        mat->cells[co.y][co.x].n_clis = 0;
    }
    mat->n_nonzeros = 0;
}

// ============================================================================
// Resende's and werneck local search functions

void update_structures(
        const rundata *run, const solution *sol, int u,
        const int *phi1, const int *phi2, const int *used,
        double *loss, double *gain, fastmat *extra, int undo){
    //
    const problem *prob = run->prob;

    int fr    = sol->assigns[u];
    double d_phi1 = -problem_assig_value(prob,phi1[u],u);
    double d_phi2 = -problem_assig_value(prob,phi2[u],u);
    assert(fr>=0 && used[fr]);
    assert(d_phi2>=d_phi1);
    assert(phi1[u]==fr);

    if(!undo){
        loss[fr] += d_phi2 - d_phi1;
    }else{
        loss[fr] -= d_phi2 - d_phi1;
        assert(loss[fr]>=-1e-6);
    }

    for(int k=0;k<prob->n_facs;k++){
        int fi = run->nearly_indexes[u][k];
        if(used[fi]) continue;
        double d_fi = -problem_assig_value(prob,fi,u);
        if(d_fi >= d_phi2) break;
        if(d_fi < d_phi1){
            if(!undo){
                gain[fi] += d_phi1 - d_fi;
                fastmat_add(extra,fi,fr,d_phi2-d_phi1);
            }else{
                gain[fi] -= d_phi1 - d_fi;
                assert(gain[fi]>=-1e-6);
                fastmat_rem(extra,fi,fr,d_phi2-d_phi1);
            }
        }else{
            if(!undo){
                fastmat_add(extra,fi,fr,d_phi2-d_fi);
            }else{
                fastmat_rem(extra,fi,fr,d_phi2-d_fi);
            }
        }
    }
}

double find_best_neighboor(
        const rundata *run,
        const double *loss, const double *gain, fastmat *extra, const int *used,
        int size_increase, int size_decrease, int *out_fins, int *out_frem){
    //
    const problem *prob = run->prob;
    //
    int best_fins = NO_MOVEMENT;
    int best_frem = NO_MOVEMENT;
    double best_delta = 0;
    // Consider just simple insertions without deletions
    if(size_increase){
        for(int fi=0;fi<prob->n_facs;fi++){
            if(used[fi]) continue;
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
        for(int fr=0;fr<prob->n_facs;fr++){
            if(!used[fr]) continue;
            double delta = prob->facility_cost[fr] - loss[fr];
            if(delta > best_delta){
                best_fins = -1;
                best_frem = fr;
                best_delta = delta;
            }
        }
    }
    // Consider swaps
    for(int i=0;i<extra->n_nonzeros;i++){
        coord co = extra->nonzeros[i];
        double extrav = extra->cells[co.y][co.x].value;
        int fi = co.y;
        int fr = co.x;
        assert(extra->cells[co.y][co.x].n_clis>0);
        assert(used[fr]);
        assert(!used[fi]);
        double delta = gain[fi] - loss[fr] + extrav - prob->facility_cost[fi] + prob->facility_cost[fr];
        if(delta > best_delta){
            best_fins = fi;
            best_frem = fr;
            best_delta = delta;
        }
    }
    // Retrieve results
    *out_fins = best_fins;
    *out_frem = best_frem;
    return best_delta;
}

int solution_resendewerneck_hill_climbing(const rundata *run, solution *sol, fastmat *zeroini_mat){
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
    assert(zeroini_mat->n_nonzeros==0);
    assert(zeroini_mat->size_y==prob->n_facs);
    assert(zeroini_mat->size_x==prob->n_facs);
    fastmat *extra = zeroini_mat;
    double *gain = safe_malloc(sizeof(double)*prob->n_facs);
    double *loss = safe_malloc(sizeof(double)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        gain[i] = 0;
        loss[i] = 0;
    }
    // Whether a facility is present or not in the solution
    int *used = safe_malloc(sizeof(int)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++) used[i] = 0;
    for(int i=0;i<sol->n_facs;i++) used[sol->facs[i]] = 1;

    // Array of affected clients
    int n_affected = prob->n_clis;
    int *affected = malloc(sizeof(int)*prob->n_clis);
    for(int i=0;i<n_affected;i++){
        affected[i] = i;
    }

    // Each movement:
    int n_moves = 0;
    int best_rem, best_ins;
    double best_delta = 0;

    int *affected_mask = malloc(sizeof(int)*prob->n_clis);

    while(1){
        best_rem = NO_MOVEMENT;
        best_ins = NO_MOVEMENT;

        // Update structures for affected clients
        for(int i=0;i<n_affected;i++){
            int u = affected[i];
            update_structures(run,sol,u,phi1,phi2,used,loss,gain,extra,0);
        }

        // Check movements that are allowed and not allowed
        int allow_size_decrease = run->local_search_rem_movement && sol->n_facs>2 &&
            (prob->size_restriction_minimum==-1 || sol->n_facs>prob->size_restriction_minimum);
        int allow_size_increase = run->local_search_add_movement &&
            (prob->size_restriction_maximum==-1 || sol->n_facs<prob->size_restriction_maximum);

        best_delta = find_best_neighboor(run,loss,gain,extra,used,
                allow_size_increase,allow_size_decrease,&best_ins,&best_rem);

        // Stop when no movement results in a better solution:
        if(best_ins==NO_MOVEMENT){
            assert(best_rem==NO_MOVEMENT);
            break;
        }

        assert(best_rem<0 || used[best_rem]);
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
            update_structures(run,sol,u,phi1,phi2,used,loss,gain,extra,1);
        }
        // Pack the list of nonzeros from the extra fastmat
        fastmat_pack(extra);

        // Perform swap:
        assert(best_ins!=NO_MOVEMENT);
        assert(best_rem!=NO_MOVEMENT);

        double old_value = sol->value;
        if(best_rem!=-1){
            solution_remove(prob,sol,best_rem,phi2,affected_mask);
            used[best_rem] = 0;
        }
        if(best_ins!=-1){
            solution_add(prob,sol,best_ins,affected_mask);
            used[best_ins] = 1;
        }
        assert(sol->value>old_value);
        #ifdef DEBUG
            double error_pred = (sol->value-old_value)-best_delta;
            if(error_pred<0) error_pred *= -1;
            assert(error_pred<1e-5 || isnan(error_pred));
        #endif

        // Update phi1 and phi2
        update_phi1_and_phi2(prob,sol,best_ins,best_rem,phi1,phi2,affected_mask);

        // Count one move:
        n_moves += 1;
    }
    assert(best_delta==0);

    free(affected_mask);

    // Clean the fastmat for reuse
    fastmat_clean(extra);

    free(affected);
    free(used);
    free(loss);
    free(gain);
    free(phi2);
    free(phi1);

    return n_moves;
}
