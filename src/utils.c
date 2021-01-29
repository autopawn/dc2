#include "utils.h"

void *safe_malloc(size_t size){
    void *ptr = malloc(size);
    if(size>0 && ptr==NULL){
        fprintf(stderr,"ERROR (on %lu B malloc, p: %p, errno: %d): %s\n",size,ptr,errno,strerror(errno));
        int currRealMem, peakRealMem, currVirtMem, peakVirtMem;
        get_memory_usage(&currRealMem,&peakRealMem,&currVirtMem,&peakVirtMem);
        fprintf(stderr,"currRealMem: %d kB\n",currRealMem);
        fprintf(stderr,"peakRealMem: %d kB\n",peakRealMem);
        fprintf(stderr,"currVirtMem: %d kB\n",currVirtMem);
        fprintf(stderr,"peakVirtmem: %d kB\n",peakVirtMem);
        exit(1);
    }
    return ptr;
}

void *safe_realloc(void *original, size_t size){
    void *ptr = realloc(original,size);
    if(size>0 && ptr==NULL){
        fprintf(stderr,"ERROR (on %lu B realloc, p: %p, errno: %d): %s\n",size,ptr,errno,strerror(errno));
        int currRealMem, peakRealMem, currVirtMem, peakVirtMem;
        get_memory_usage(&currRealMem,&peakRealMem,&currVirtMem,&peakVirtMem);
        fprintf(stderr,"currRealMem: %d kB\n",currRealMem);
        fprintf(stderr,"peakRealMem: %d kB\n",peakRealMem);
        fprintf(stderr,"currVirtMem: %d kB\n",currVirtMem);
        fprintf(stderr,"peakVirtmem: %d kB\n",peakVirtMem);
        exit(1);
    }
    return ptr;
}

// Thanks to user Thomas Mueller: https://stackoverflow.com/a/12996028
uint hash_int(uint x){
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

int elem_in_sorted(int *array, int len, int val){
    int a = 0;
    int b = len;
    while(a<b){
        int c = (a+b)/2;
        if(array[c]==val) return 1;
        if(val<array[c]){
            b = c;
        }else{
            a = c+1;
        }
    }
    return 0;
}

int diff_sorted(int *arr1, int len1, int *arr2, int len2){
    int diff = 0;
    int i1 = 0;
    int i2 = 0;
    while(i1<len1 && i2<len2){
        int v1 = arr1[i1];
        int v2 = arr2[i2];
        if(v1==v2){
            i1 += 1;
            i2 += 1;
        }else if(v1<v2){
            diff += 1;
            i1 += 1;
        }else{
            diff += 1;
            i2 += 1;
        }
    }
    diff += (i1<i2)? i2-i1 : i1-i2;
    return diff;
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// Thanks to user Anti Earth: https://stackoverflow.com/a/47531152
/*
 * Measures the current (and peak) resident and virtual memories
 * usage of your linux C process, in kB
 */
void get_memory_usage(int* currRealMem, int* peakRealMem, int* currVirtMem, int* peakVirtMem){
    // stores each word in status file
    static char memory_buffer[1024];
    memory_buffer[0] = '\0';

    // Initialize values to -1 if couldn't be read.
    if(currRealMem) *currRealMem = -1;
    if(peakRealMem) *peakRealMem = -1;
    if(currVirtMem) *currVirtMem = -1;
    if(peakVirtMem) *peakVirtMem = -1;

    // linux file contains this-process info
    FILE* file = fopen("/proc/self/status", "r");
    if(file==NULL){
        fprintf(stderr,"WARNING: couldn't read /proc/self/status to retrieve memory usage.\n");
        return;
    }

    // read the entire file
    while (fscanf(file, " %1023s", memory_buffer) == 1) {

        if (currRealMem && strcmp(memory_buffer, "VmRSS:") == 0) {
            if (fscanf(file, " %d", currRealMem)!=1) *currRealMem = -1;
        }
        if (peakRealMem && strcmp(memory_buffer, "VmHWM:") == 0) {
            if (fscanf(file, " %d", peakRealMem)!=1) *peakRealMem = -1;
        }
        if (currVirtMem && strcmp(memory_buffer, "VmSize:") == 0) {
            if (fscanf(file, " %d", currVirtMem)!=1) *currVirtMem = -1;
        }
        if (peakVirtMem && strcmp(memory_buffer, "VmPeak:") == 0) {
            if (fscanf(file, " %d", peakVirtMem)!=1) *peakVirtMem = -1;
        }
    }
    fclose(file);
}
