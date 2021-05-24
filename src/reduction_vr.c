#include "reduction.h"

typedef long long int lint;

// The dissimilitude between solutions whit ids, id1 and id2.
typedef struct {
    double dissim;
    int id1, id2;
} dissimpair;


// Compare dissimpairs by dissimilitude (<0 means a has less dissimilitude)
int dissimpair_cmp(dissimpair a, dissimpair b){
    double delta = a.dissim - b.dissim;
    if(delta!=0) return delta;
    delta = a.id1-b.id1;
    if(delta!=0) return delta;
    delta = a.id2-b.id2;
    return (delta==0)? 0 : ((delta<0)? -1 : +1);
}

// A heap of dissimpairs
typedef struct {
    dissimpair *elems;
    lint len;
    lint size;
} pairheap;

// Initializes a heap of dissimpairs.
pairheap *pairheap_init(lint size){
    pairheap *heap = safe_malloc(sizeof(pairheap));
    heap->size = size;
    heap->len = 0;
    heap->elems = safe_malloc(sizeof(dissimpair)*heap->size);
    return heap;
}

// Free the memory required by a heap of dissimpairs.
void pairheap_free(pairheap *heap){
    free(heap->elems);
    free(heap);
}

// Gets the dissimpair of smaller dissimilitude of the heap, removing it.
dissimpair pairheap_poll(pairheap *heap){
    assert(heap->len>0);
    dissimpair ret = heap->elems[0];
    heap->elems[0] = heap->elems[heap->len-1];
    heap->len -= 1;
    // Heapify down:
    lint i = 0;
    lint c = 2*i+1;
    while(c<heap->len){
        if(c+1<heap->len && dissimpair_cmp(heap->elems[c+1],heap->elems[c])<0) c = c+1;
        if(dissimpair_cmp(heap->elems[i],heap->elems[c])<0) break;
        dissimpair aux = heap->elems[i];
        heap->elems[i] = heap->elems[c];
        heap->elems[c] = aux;
        i = c;
        c = 2*i+1;
    }
    return ret;
}

// Adds a dissimpair to the heap.
void pairheap_add(pairheap *heap, dissimpair val){
    assert(heap->len<heap->size);
    heap->len += 1;
    heap->elems[heap->len-1] = val;
    // Heapify up:
    lint i = heap->len-1;
    lint p = (i-1)/2;
    while(i>0 && dissimpair_cmp(heap->elems[i],heap->elems[p])<0){
        dissimpair aux = heap->elems[i];
        heap->elems[i] = heap->elems[p];
        heap->elems[p] = aux;
        i = p;
        p = (i-1)/2;
    }
}

typedef struct {
    // Parameters
    int thread_id;
    int vision_range;
    soldismode soldis;
    facdismode facdis;
    // Concurrent axes to heap
    pthread_mutex_t *heap_mutex;
    pairheap *heap;
    // Solutions (read only)
    const rundata *run;
    int n_sols;
    const solution **sols;
    // Fast way of getting previous and next solutions
    const int *prev_sols;
    const int *next_sols;
    // Semaphores
    sem_t *thread_sem;   // The main thread finished.
    sem_t *complete_sem; // Inform that thread has finished.
    int *terminated;     // For the main thread to inform that the job is done.
} reductionvr_thread_args;

void *reductionvr_thread_execution(void *arg){
    reductionvr_thread_args *args = (reductionvr_thread_args *) arg;

    { // Help building initial set of dissimilitude pairs
        dissimpair *pairs = safe_malloc(sizeof(dissimpair)*args->vision_range);
        int n_pairs = 0;
        for(int i=args->thread_id;i<args->n_sols;i+=args->run->n_threads){
            for(int j=1;j<=args->vision_range;j++){
                if(i+j>=args->n_sols) break;
                pairs[n_pairs].id1 = i;
                pairs[n_pairs].id2 = i+j;
                pairs[n_pairs].dissim = solution_dissimilitude(
                    args->run,args->sols[i],args->sols[i+j],
                    args->soldis,args->facdis);
                n_pairs += 1;
            }
            // Add pairs to the heap
            pthread_mutex_lock(args->heap_mutex);
                for(int k=0;k<n_pairs;k++) pairheap_add(args->heap,pairs[k]);
            pthread_mutex_unlock(args->heap_mutex);
            // Delete pairs
            n_pairs = 0;
        }
        free(pairs);
    }

    { // Wake when required to compute new dissimilitudes
        dissimpair *pair_buffer = safe_malloc(sizeof(dissimpair)
            *(1+args->vision_range/args->run->n_threads));
        int pair_buffer_len = 0;
        // Each time a solution is deleted
        while(1){
            sem_post(args->complete_sem);
            // ---@> Main thread works here.
            sem_wait(args->thread_sem);
            assert(args->prev_sols!=NULL);
            assert(args->next_sols!=NULL);
            // If main thread says the job is done, break
            if(*args->terminated) break;
            // Create new pairs
            for(int i=args->thread_id;i<args->vision_range;i+=args->run->n_threads){
                int pair_a = args->prev_sols[args->vision_range-1-i];
                int pair_b = args->next_sols[i];
                if(pair_a!=-1 && pair_b!=-1){
                    // Create the replace dissimpair:
                    dissimpair pair;
                    pair.id1 = pair_a;
                    pair.id2 = pair_b;
                    pair.dissim = solution_dissimilitude(args->run,
                        args->sols[pair_a],args->sols[pair_b],
                        args->soldis,args->facdis);
                    // Add to pair buffer:
                    pair_buffer[pair_buffer_len] = pair;
                    pair_buffer_len += 1;
                    // Add pairs in buffer to to the heap if mutex is free
                    if(pthread_mutex_trylock(args->heap_mutex)==0){
                        for(int k=0;k<pair_buffer_len;k++){
                            pairheap_add(args->heap,pair_buffer[k]);
                        }
                        pthread_mutex_unlock(args->heap_mutex);
                        pair_buffer_len = 0;
                    }
                }
            }
            // If pairs remain, wait for mutex
            if(pair_buffer_len>0){
                pthread_mutex_lock(args->heap_mutex);
                for(int i=0;i<pair_buffer_len;i++){
                    pairheap_add(args->heap,pair_buffer[i]);
                }
                pthread_mutex_unlock(args->heap_mutex);
                pair_buffer_len = 0;
            }
        }
        free(pair_buffer);
    }
    return NULL;
}

void reduction_vr_heuristic(const rundata *run, solution **sols, int *n_sols,
        int n_target, soldismode soldis, facdismode facdis, int vision_range){
    assert(vision_range>0);

    // Ensure that the vision_range isn't larger than the number of solutions.
    if(vision_range>*n_sols) vision_range = *n_sols;
    // Sort solutions in decreasing order
    qsort(sols,*n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    // Return if there is no need of reduction.
    if(*n_sols<=n_target) return;
    // To know if a solution has been discarded:
    int *discarted = safe_malloc((*n_sols)*sizeof(int));
    for(int i=0;i<*n_sols;i++) discarted[i] = 0;

    // The solutions previous and next to the last deleted solution
    // (to communicate with the threads).
    int *prev_sols = safe_malloc(sizeof(int)*vision_range);
    int *next_sols = safe_malloc(sizeof(int)*vision_range);

    // Heap of dissimilitude pairs
    pairheap *heap = pairheap_init(2*(*n_sols)*vision_range);
    pthread_mutex_t heap_mutex;
    pthread_mutex_init(&heap_mutex,NULL);

    // Create threads:
    pthread_t *threads = safe_malloc(sizeof(pthread_t)*run->n_threads);
    sem_t **t_sems = safe_malloc(sizeof(sem_t *)*run->n_threads);
    sem_t **c_sems = safe_malloc(sizeof(sem_t *)*run->n_threads);
    reductionvr_thread_args *targs = safe_malloc(sizeof(reductionvr_thread_args)*run->n_threads);
    int terminated = 0; // NOTE: Only access before threads are waken.
    //
    for(int i=0;i<run->n_threads;i++){
        // Init semaphores
        t_sems[i] = dc_semaphore_init();
        c_sems[i] = dc_semaphore_init();

        // Parameters
        targs[i].thread_id = i;
        targs[i].vision_range = vision_range;
        targs[i].soldis = soldis;
        targs[i].facdis = facdis;
        // Concurrent axes to heap
        targs[i].heap_mutex = &heap_mutex;
        targs[i].heap = heap;
        // Solutions (read only)
        targs[i].run = run;
        targs[i].n_sols = *n_sols;
        targs[i].sols = (const solution **) sols;
        // Fast way of getting previous and next solutions
        targs[i].prev_sols = prev_sols;
        targs[i].next_sols = next_sols;
        // Semaphores
        targs[i].thread_sem = t_sems[i];
        targs[i].complete_sem = c_sems[i];
        targs[i].terminated = &terminated;

        // Initialize thread
        int rc = pthread_create(&threads[i],NULL,reductionvr_thread_execution,&targs[i]);
        if(rc){
            fprintf(stderr,"Error %d on thread creation\n",rc);
            exit(1);
        }
    }

    // Double linked list to know which solution comes before and after it
    int *nexts = safe_malloc((*n_sols)*sizeof(int));
    int *prevs = safe_malloc((*n_sols)*sizeof(int));
    for(int i=0;i<*n_sols;i++){
        prevs[i] = i-1;
        nexts[i] = i+1;
    }
    prevs[0] = -1;
    nexts[*n_sols-1] = -1;

    // Wait for all threads to terminate the last job
    for(int i=0;i<run->n_threads;i++){
        sem_wait(c_sems[i]);
    }
    // ---@> At this point, all initial dissimilitude pairs are complete.

    // Eliminate as many solutions as required:
    int n_eliminate = *n_sols-n_target;
    int elims = 0;
    while(elims<n_eliminate){
        // Eliminate worst solution of most similar pair
        if(heap->len==0) break;
        dissimpair pair = pairheap_poll(heap);
        if(!discarted[pair.id1] && !discarted[pair.id2]){
            // Delete the second solution on the pair.
            int to_delete = pair.id2;
            discarted[to_delete] = 1;
            solution_free(sols[to_delete]);
            elims += 1;
            // Update double linked list:
            if(nexts[to_delete]!=-1) prevs[nexts[to_delete]] = prevs[to_delete];
            if(prevs[to_delete]!=-1) nexts[prevs[to_delete]] = nexts[to_delete];
            // Add new pairs to replace those that will be deleted on the destroyed solution.
            int iter;
            // Get solutions after
            iter = to_delete;
            for(int i=0;i<vision_range;i++){
                if(nexts[iter]==-1){
                    next_sols[i] = -1;
                }else{
                    iter = nexts[iter];
                    next_sols[i] = iter;
                }
            }
            // Get solutions before
            iter = to_delete;
            for(int i=0;i<vision_range;i++){
                if(prevs[iter]==-1){
                    prev_sols[i] = -1;
                }else{
                    iter = prevs[iter];
                    prev_sols[i] = iter;
                }
            }
            // Wake threads to create new pairs
            for(int i=0;i<run->n_threads;i++){
                sem_post(t_sems[i]);
            }
            // Wait for all threads to terminate their job
            for(int i=0;i<run->n_threads;i++){
                sem_wait(c_sems[i]);
            }
        }
    }

    // Wake threads for termination
    terminated = 1;
    for(int i=0;i<run->n_threads;i++){
        sem_post(t_sems[i]);
    }

    // Join threads
    for(int i=0;i<run->n_threads;i++){
        pthread_join(threads[i],NULL);
    }
    free(threads);
    free(targs);

    // Destroy semaphores
    for(int i=0;i<run->n_threads;i++){
        dc_semaphore_free(t_sems[i]);
        dc_semaphore_free(c_sems[i]);
    }
    free(t_sems);
    free(c_sems);

    // Free prev and next sols arrays
    free(prev_sols);
    free(next_sols);
    // Free all the pairs:
    pairheap_free(heap);
    // Set output final array:
    int new_nsols=0;
    for(int i=0;i<*n_sols;i++){
        if(discarted[i]==0){
            sols[new_nsols] = sols[i];
            new_nsols += 1;
        }
    }
    *n_sols = new_nsols;
    // Free arrays
    free(discarted);
    free(nexts);
    free(prevs);
}
