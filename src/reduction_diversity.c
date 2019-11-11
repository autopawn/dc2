#include "reduction.h"

typedef struct {
    int thread_id;
    int n_target;
    int *centroids;
    int *nearest_cluster;
    double *nearest_dist;
    const problem *prob;
    int n_sols;
    solution **sols;
    soldismode soldis;
    facdismode facdis;
    sem_t *thread_sem;   // The main thread picked the next centroid.
    sem_t *complete_sem; // Inform that thread has terminated updating its part of nearest_cluster.
} reductiondiv_thread_args;

void *reductiondiv_thread_execution(void *arg){
    reductiondiv_thread_args *args = (reductiondiv_thread_args *) arg;
    for(int t=0;t<args->n_target;t++){
        sem_wait(args->thread_sem);
        int centroid = args->centroids[t];
        // Help updating the dissimilitudes to each cluster
        for(int r=args->thread_id;r<args->n_sols;r+=args->prob->n_threads){
            double disim = solution_dissimilitude(args->prob,
                args->sols[r],args->sols[centroid],args->soldis,args->facdis);
            if(disim<args->nearest_dist[r]){
                args->nearest_cluster[r] = t;
                args->nearest_dist[r] = disim;
            }
        }
        sem_post(args->complete_sem);
    }
    return NULL;
}

void reduction_diversity_starting(const problem *prob, solution **sols, int *n_sols,
        int n_target, soldismode soldis, facdismode facdis, int bests_of_clusters){
    // Sort solutions in decreasing order
    qsort(sols,*n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    if(*n_sols<=n_target) return;
    // Indexes of selected centroids
    int *centroids = safe_malloc(sizeof(int)*n_target);
    int *is_centroid = safe_malloc(sizeof(int)*(*n_sols));
    centroids[0] = 0; 
    int n_centroids = 1;
    // Nearest cluster index and distance to it
    int *nearest_cluster = safe_malloc(sizeof(int)*(*n_sols));
    double *nearest_dist = safe_malloc(sizeof(double)*(*n_sols));
    // Creation of the first cluster including all the solutions.
    for(int i=0;i<*n_sols;i++){
        nearest_cluster[i] = -1;
        nearest_dist[i] = INFINITY;
        is_centroid[i] = 0;
    }
    is_centroid[0] = 1;
    // Create threads
    pthread_t *threads = safe_malloc(sizeof(pthread_t)*prob->n_threads);
    sem_t *t_sems = safe_malloc(sizeof(sem_t)*prob->n_threads);
    sem_t *c_sems = safe_malloc(sizeof(sem_t)*prob->n_threads);
    reductiondiv_thread_args *targs = safe_malloc(sizeof(reductiondiv_thread_args)*prob->n_threads);
    for(int i=0;i<prob->n_threads;i++){
        // Init semaphores
        if(sem_init(&t_sems[i],0,0)==-1 || sem_init(&c_sems[i],0,0)==-1){
            fprintf(stderr,"ERROR: sem_init failed!\n");
            exit(1);
        }
        targs[i].thread_id = i;
        targs[i].n_target = n_target;
        targs[i].centroids = centroids;
        targs[i].nearest_cluster = nearest_cluster;
        targs[i].nearest_dist = nearest_dist;
        targs[i].prob = prob;
        targs[i].n_sols = *n_sols;
        targs[i].sols = sols;
        targs[i].soldis = soldis;
        targs[i].facdis = facdis;
        targs[i].thread_sem = &t_sems[i];
        targs[i].complete_sem = &c_sems[i];
        int rc = pthread_create(&threads[i],NULL,reductiondiv_thread_execution,&targs[i]);
        if(rc){
            fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
            exit(1);
        }
    }
    for(int t=0;t<n_target;t++){
        // Allow threads to update the distances again
        for(int i=0;i<prob->n_threads;i++) sem_post(&t_sems[i]);
        // Wait for all threads to have updated the centroids
        for(int i=0;i<prob->n_threads;i++) sem_wait(&c_sems[i]);
        
        if(t==n_target-1) break;

        // Find the new centroid
        int farthest = -1;
        double fardist = -INFINITY;
        for(int i=1;i<*n_sols;i++){
            if(!is_centroid[i] && nearest_dist[i] > fardist){
                farthest = i;
                fardist = nearest_dist[farthest];
            }
        }
        assert(farthest!=-1);

        // Add the new centroid
        centroids[n_centroids] = farthest;
        assert(!is_centroid[farthest]);
        is_centroid[farthest] = 1;
        n_centroids += 1;
    }
    assert(n_centroids==n_target);

    /* Free thread memory and close semaphores. */
    for(int i=0;i<prob->n_threads;i++){ // Join threads
        pthread_join(threads[i],NULL);
    }
    free(targs);
    for(int i=0;i<prob->n_threads;i++){ // Close semaphores
        if(sem_destroy(&c_sems[i])==-1 || sem_destroy(&t_sems[i])==-1){
            fprintf(stderr,"ERROR: sem_destroy failed!\n");
            exit(1);
        }
    }
    free(c_sems);
    free(t_sems);
    free(threads);

    /* Select final solutions */
    // Ensure that the centroid is part of its cluster (in case of distance 0 ties)
    for(int i=0;i<n_centroids;i++){
        assert(is_centroid[centroids[i]]);
        nearest_cluster[centroids[i]] = i;
    }
    // Replace the centroids with the best of each cluster
    if(bests_of_clusters){
        // Reset all centroids:
        for(int i=0;i<n_centroids;i++) centroids[i] = -1;
        // For each solution check if is the best of its cluster
        for(int i=0;i<*n_sols;i++){
            // Check if it is new centroid (first found on the given cluster).
            if(centroids[nearest_cluster[i]]==-1){
                centroids[nearest_cluster[i]] = i;
                is_centroid[i] = 1;
            }else{
                is_centroid[i] = 0;
            }
        }
    }
    // Select all solutions that are centroids
    int n_final = 0;
    for(int i=0;i<*n_sols;i++){
        if(is_centroid[i]){
            sols[n_final] = sols[i]; 
            n_final += 1;
        }else{
            solution_free(sols[i]);
        }
    }
    assert(n_final==n_target);
    *n_sols = n_target;

    // Free remaining arrays
    free(nearest_dist);
    free(nearest_cluster);
    free(is_centroid);
    free(centroids);
}