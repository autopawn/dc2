#include "fastmat.h"
#include "utils.h"

#include <string.h>
#include <stdio.h>

// Initialize a fastmat
fastmat *fastmat_init(){
    fastmat *mat = malloc(sizeof(fastmat));
    mat->tab = coordcelltable_init();
    return mat;
}

// Free fastmat from memory
void fastmat_free(fastmat *mat){
    coordcelltable_free(mat->tab);
    free(mat);
}

// Adds a value on the given position in the fastmatrix
void fastmat_add(fastmat *mat, int y, int x, double v){
    // Add cell to the fastmatrix if not already there
    if(!coordcelltable_has(mat->tab,(coord){.x=x,.y=y})){
        coordcelltable_add(mat->tab,(coord){.x=x,.y=y},(cell){.value=0,.n_clis=0});
    }
    // Modify the current cell, add the value
    cell c = coordcelltable_get(mat->tab,(coord){.x=x,.y=y});
    c.value += v;
    c.n_clis += 1;
    // Update in the coordcelltable
    coordcelltable_add(mat->tab,(coord){.x=x,.y=y},c);
}

// Intended for reverting an addition on the given position in a fastmatrix
void fastmat_rem(fastmat *mat, int y, int x, double v){
    assert(coordcelltable_has(mat->tab,(coord){.x=x,.y=y}));
    // Substract and addition to the cell
    cell c = coordcelltable_get(mat->tab,(coord){.x=x,.y=y});
    c.value -= v;
    c.n_clis -= 1;
    if(c.n_clis==0){
        // Delete cell
        coordcelltable_pop(mat->tab,(coord){.x=x,.y=y});
    }else{
        // Update cell
        coordcelltable_add(mat->tab,(coord){.x=x,.y=y},c);
    }
}

// Resets a fastmat (so that it only contains zeroes)
void fastmat_clean(fastmat *mat){
    coordcelltable_free(mat->tab);
    mat->tab = coordcelltable_init();
}

double fastmat_get(fastmat *mat, int x, int y){
    if(coordcelltable_has(mat->tab,(coord){.x=x,.y=y})){
        return coordcelltable_get(mat->tab,(coord){.x=x,.y=y}).value;
    }
    return 0.0;
}






/*


// Initialize a fastmat
fastmat *fastmat_init(){
    fastmat *mat = safe_malloc(sizeof(fastmat));

    // Initialize hashtable
    memset(mat->slots,0,sizeof(mat->slots));
    memset(mat->n_slots,0,sizeof(mat->n_slots));

    // Initialize dynamic array
    mat->s_nonzeros = INITIAL_NONZEROARRAY_SIZE;
    mat->n_nonzeros = 0;
    mat->nonzeros = safe_malloc(sizeof(cell)*mat->n_nonzeros);

    return mat;
}

// Find nonzero cell index
int find_idx(fastmat *mat, int y, int x){
    int slot = HASHPOS(x,y);
    for(int i=0;i<mat->n_slots[slot];i++){
        coord c = mat->slots[slot][i];
        if(c.x == x && c.y == y) return c.i;
    }
    return -1;
}

// Adds a value on the given position in the fastmatrix
void fastmat_add(fastmat *mat, int y, int x, double v){
    int idx = find_idx(mat,y,x);
    // If nonzero cell doesn't exist yet, create it
    if(idx==-1){
        // Extend dynamic array
        if(mat->n_nonzeros==mat->s_nonzeros){
            mat->s_nonzeros *= 2;
            safe_realloc(mat->nonzeros,sizeof(cell)*mat->s_nonzeros);
        }

        // New nonzero cell (also set idx)
        idx = mat->n_nonzeros;
        mat->n_nonzeros += 1;
        mat->nonzeros[idx].value  = 0;
        mat->nonzeros[idx].n_clis = 0;

        // New coord in hashtable
        int slot = HASHPOS(x,y);
        int j = mat->n_slots[slot];
        mat->n_slots[slot] += 1;
        safe_realloc(mat->slots[slot],sizeof(coord)*mat->n_slots[slot]);
        mat->slots[slot][j].x = x;
        mat->slots[slot][j].y = y;
        mat->slots[slot][j].i = idx;
    }
    // Add to nonzero cell
    assert(idx>=0);
    assert(idx<mat->n_nonzeros);
    assert(v>=-1e-6);
    mat->nonzeros[idx].value += v;
    mat->nonzeros[idx].n_clis += 1;
}

// Intended for reverting an addition on the given position in a fastmatrix
void fastmat_rem(fastmat *mat, int y, int x, double v){
    int idx = find_idx(mat,y,x);
    assert(idx>=0);
    mat->nonzeros[idx].value -= v;
    mat->nonzeros[idx].n_clis -= 1;
    // Destroy nonzero if reaches 0 aportations
    if(mat->nonzeros[idx].n_clis==0){
        #ifdef DEBUG
            assert(mat->nonzeros[idx].value<=0.000001);
        #endif

        // == Delete nonzero, replace it with last nonzero
        int n_nonz = mat->n_nonzeros;
        mat->nonzeros[idx].n_clis = mat->nonzeros[n_nonz-1].n_clis;
        mat->nonzeros[idx].value  = mat->nonzeros[n_nonz-1].value;
        mat->n_nonzeros -= 1;

        // == Update coord of last nonzero
        int slot = HASHPOS(x,y);


        // == Delete the coord in hashtable
        int slot = HASHPOS(x,y);
        int n_slots = mat->n_slots[slot];
        // Find and replace coord with last coord in this slot
        for(int i=0;i<n_slots;i++){
            coord *c = &mat->slots[slot][i];
            if(c->x == x && c->y == y){
                c->x = mat->slots[slot][n_slots-1].x;
                c->y = mat->slots[slot][n_slots-1].y;
                break;
            }
        }
        // Make slot smaller
        mat->n_slots[slot] -= 1;
        safe_realloc(mat->slots[slot],sizeof(coord)*mat->n_slots[slot]);
    }
}

// Free fastmat from memory
void fastmat_free(fastmat *mat);

// Resets a fastmat (so that it only contains zeroes)
void fastmat_clean(fastmat *mat);















// Initializes fastmat
fastmat *fastmat_init(int size_y, int size_x){
    fastmat *mat = safe_malloc(sizeof(fastmat));
    mat->cells = safe_malloc(sizeof(cell*)*size_y);
    for(int y=0;y<size_y;y++){
        mat->cells[y] = safe_malloc(sizeof(cell)*size_x);
        for(int x=0;x<size_x;x++){
            mat->cells[y][x].value  = 0;
            mat->cells[y][x].n_clis = 0;
            mat->cells[y][x].nonzero_index = -1;
        }
    }
    mat->size_y = size_y;
    mat->size_x = size_x;
    mat->s_nonzeros = size_y*size_x;
    mat->nonzeros = safe_malloc(sizeof(coord)*mat->s_nonzeros);
    mat->n_nonzeros = 0;
    return mat;
}

// Adds a value on the given position in the fastmatrix
void fastmat_add(fastmat *mat, int y, int x, double v){
    // Find current cell
    assert(y>=0 && x>=0 && y<mat->size_y && x<mat->size_x);
    assert(v>=-1e-6);
    mat->cells[y][x].value += v;
    mat->cells[y][x].n_clis += 1;
    if(mat->cells[y][x].n_clis==1){ // Add nonzero to the list
        assert(mat->n_nonzeros<mat->s_nonzeros);
        coord *co = &mat->nonzeros[mat->n_nonzeros];
        co->y = y;
        co->x = x;
        mat->cells[y][x].nonzero_index = mat->n_nonzeros;
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

        // Delete this entry from the nonzeros array, swap with the last nonzero
        int c_index = mat->cells[y][x].nonzero_index;

        mat->nonzeros[c_index] = mat->nonzeros[mat->n_nonzeros-1];
        int co_y = mat->nonzeros[c_index].y;
        int co_x = mat->nonzeros[c_index].x;
        mat->cells[co_y][co_x].nonzero_index = c_index;

        mat->cells[y][x].nonzero_index = -1;

        mat->n_nonzeros--;
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

*/
