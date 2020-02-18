#include "reduction.h"

typedef long long int weight;

// Updates the probability weight at position k in the given segmentree.
void segtree_update(weight *segtree, int len, int k, weight val){
    int p = len + k;
    weight delta = val - segtree[p];
    while(p>0){
        segtree[p] += delta;
        p /= 2;
    }
}

// Picks a node at random from the segment tree, retrieves it.
int segtree_pick(weight *segtree, int len){
    weight total = segtree[1];
    assert(total>0);
    weight r = (weight)(((double)rand()/RAND_MAX)*total);
    int p = 1;
    while(p<len){
        if(r>=segtree[2*p]){
            r -= segtree[2*p];
            p = 2*p+1;
        }else{
            p = 2*p;
        }
    }
    return p-len;
}

void reduction_random_rank(const rundata *run, solution **sols, int *n_sols, int n_target, int elitist){
    if(*n_sols<=n_target) return;
    // Sort solutions in decreasing order
    qsort(sols,*n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    // Solutions that were picked
    int *picked = safe_malloc(sizeof(int)*(*n_sols));
    for(int i=0;i<(*n_sols);i++) picked[i] = 0;
    // Probability weight segment tree
    weight *pws = safe_malloc(sizeof(weight)*(*n_sols)*2);
    for(int i=0;i<(*n_sols)*2;i++){
        pws[i] = 0;
    }
    // The first solution is always picked if elitist (unless n_target==0)
    picked[0] = (elitist && n_target>0);
    int n_picked = picked[0];
    // Initial probabilities
    for(int i=elitist;i<*n_sols;i++){
        weight pw = 100000000/(i+1);
        if(pw<1) pw = 1;
        segtree_update(pws,*n_sols,i,pw);
    }
    // Start picking solutions
    while(n_picked<n_target){
        // Get a solution at random from the segment tree
        int k = segtree_pick(pws,*n_sols);
        picked[k] = 1;
        segtree_update(pws,*n_sols,k,0);
        n_picked++;
    }
    // Delete not picked solutions
    int n_final = 0;
    for(int i=0;i<*n_sols;i++){
        if(picked[i]){
            sols[n_final] = sols[i];
            n_final++;
        }else{
            solution_free(sols[i]);
        }
    }
    *n_sols = n_final;
    // Free memory
    free(pws);
    free(picked);
}