#include "shuffle.h"

#include <assert.h>

#define RANXOSHI256_IMPLEMENTATION
#include "lib/ranxoshi256.h"


shuffler *shuffler_init(int len){
    // Allocate memory and initialize numbers
    shuffler *shu = safe_malloc(sizeof(shuffler));
    shu->len = len;
    shu->nums = safe_malloc(sizeof(unsigned int)*shu->len);
    for(int i=0;i<len;i++) shu->nums[i] = i;
    // Initialize random number generator
    unsigned char seed[32];
    for(int i=0;i<32;i++) seed[i] = rand()%256;
    ranxoshi256Seed(&shu->rngen,seed);
    // First shuffler, shuffler everything
    shu->n_retrieved = shu->len;
    shuffler_reshuffle(shu);
    //
    return shu;
}

unsigned int shuffler_next(shuffler *shu){
    assert(shu->n_retrieved<shu->len);
    unsigned int v = shu->nums[shu->n_retrieved];
    shu->n_retrieved += 1;
    return v;
}

void shuffler_reshuffle(shuffler *shu){
    for(int i=0;i<shu->n_retrieved;i++){
        int repl = i+(ranxoshi256Next(&shu->rngen)%(shu->len-i));
        int aux = shu->nums[i];
        shu->nums[i] = shu->nums[repl];
        shu->nums[repl] = aux;
    }
    shu->n_retrieved = 0;
}

void shuffler_free(shuffler *shu){
    free(shu->nums);
    free(shu);
}
