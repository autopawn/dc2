#include "utils.h"

void *safe_malloc(size_t size){
    assert(errno==0);
    void *ptr = malloc(size);
    if((size>0 && ptr==NULL) || errno!=0){
        fprintf(stderr,"ERROR (on malloc): %s\n",strerror(errno));
        exit(1);
    }
    return ptr;
}

uint hash_int(uint x){
    // Thanks to https://stackoverflow.com/a/12996028
    x = ((x >> 16)^x)*0x45d9f3b;
    x = ((x >> 16)^x)*0x45d9f3b;
    x = (x >> 16)^x;
    return x;
}

void add_to_sorted(int *array, int *len, int val){
    int place = *len;
    while(place>0){
        if(array[place-1]<=val) break;
        array[place] = array[place-1];
        place--;
    }
    array[place] = val;
    *len += 1;
}

void rem_of_sorted(int *array, int *len, int val){
    int place=-1;
    for(int i=0;i<*len;i++){
        if(array[i]==val){
            place = i;
            break;
        }
    }
    assert(place!=-1); // Not in array.
    for(int i=place;i<*len-1;i++){
        array[i] = array[i+1];
    }
    *len -= 1;
}

sem_t *dc_semaphore_init(){
    sem_t *sem;
    assert(errno==0);
    #ifdef NAMED_SEMAPHORES
        // id changes each time this function is called
        static int sem_id = 1000;
        char namebuffer[200];
        sprintf(namebuffer,"/dc2_sem_%d",sem_id++);
        sem = sem_open(namebuffer,O_CREAT | O_EXCL,0644,0);
        if(sem == SEM_FAILED){
            perror("ERROR: sem_open failed!");
            exit(1);
        }
        if(sem_unlink(namebuffer)!=0){
            perror("ERROR: sem_unlink failed!");
            exit(1);
        }
    #else
        sem = safe_malloc(sizeof(sem_t)*1);
        if(sem_init(sem,0,0)==-1){
            perror("ERROR: sem_init failed!");
            exit(1);
        }
    #endif
    assert(errno==0);
    return sem;
}

void dc_semaphore_free(sem_t *sem){
    assert(errno==0);
    #ifdef NAMED_SEMAPHORES
        if(sem_close(sem)==-1){
            perror("ERROR: sem_close failed!");
            exit(1);
        }
        assert(errno==0);
    #else
        if(sem_destroy(sem)==-1){
            perror("ERROR: sem_destroy failed!");
            exit(1);
        }
        free(sem);
    #endif
    assert(errno==0);
}