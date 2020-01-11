#include "load.h"

problem *load_simple_format(FILE *fp){
    // Read filename in header:
    char buffer[400];
    if(fscanf(fp,"FILE: %s",buffer)!=1){
        fprintf(stderr,"ERROR: couldn't read FILE!\n");
        exit(1);
    }

    // Read the number of facilities:
    int n_facs;
    if(fscanf(fp,"%d",&n_facs)!=1){
        fprintf(stderr,"ERROR: number of facilities expected!\n");
        exit(1);
    }

    // Read the number of clients:
    int n_clis;
    if(fscanf(fp,"%d",&n_clis)!=1){
        fprintf(stderr,"ERROR: number of clients expected!\n");
        exit(1);
    }

    // Allocate problem
    problem *prob = problem_init(n_facs,n_clis);

    // Third argument is the size restriction.
    int trd_num;
    if(fscanf(fp,"%d",&trd_num)!=1){
        fprintf(stderr,"ERROR: size restriction expected!\n");
        exit(1);
    }
    if(trd_num==0) trd_num=-1;
    prob->size_restriction_maximum = trd_num;
    prob->size_restriction_minimum = trd_num;

    // For each facility
    for(int i=0;i<prob->n_facs;i++){

        // Read facility index
        int facility_index;
        if(fscanf(fp,"%d",&facility_index)!=1){
            fprintf(stderr,"ERROR: facility index expected!\n");
        }
        assert(i+1==facility_index);

        // Read facility cost
        if(fscanf(fp,"%lf",&prob->facility_cost[i])!=1){
            fprintf(stderr,"ERROR: facility %d cost expected!\n",i);
            exit(1);
        }

        // Read each distance
        for(int j=0;j<prob->n_clis;j++){
            if(fscanf(fp,"%lf",&prob->distance[i][j])!=1){
                fprintf(stderr,"ERROR: distance from facility %d to client %d expected!\n",i,j);
                exit(1);
            }
        }
    }

    // Unsetted values (for SPLP)
    for(int j=0;j<prob->n_clis;j++){
        prob->client_weight[j] = 1;
    }

    prob->transport_cost = 1;
    prob->client_gain = 0;
    prob->unassigned_cost = INFINITY;
    prob->filter = BETTER_THAN_ALL_PARENTS;
    prob->lower_bound = -INFINITY;
    prob->branch_and_bound = 1;

    //
    return prob;
}

problem *load_dc_format(FILE *fp, int dc_version){
    // For now, just this one version
    assert(dc_version==1);

    // Read line in header
    char buffer[400];
    if(fscanf(fp,"%s",buffer)!=1){
        fprintf(stderr,"ERROR: couldn't read header!\n");
        exit(1);
    }

    // Parameters
    double transport_cost;
    double client_gain;
    double unassigned_cost;
    int min_size,max_size;
    int n_facs,n_clis;

    // Scan parameters
    if(fscanf(fp,"%lf",&transport_cost)!=1){
        fprintf(stderr,"ERROR: transport cost expected!\n");
        exit(1);
    }
    if(fscanf(fp,"%lf",&client_gain)!=1){
        fprintf(stderr,"ERROR: client gain expected!\n");
        exit(1);
    }
    if(fscanf(fp,"%lf",&unassigned_cost)!=1){
        fprintf(stderr,"ERROR: unassigned cost expected!\n");
        exit(1);
    }
    if(fscanf(fp,"%d",&min_size)!=1){
        fprintf(stderr,"ERROR: minimum size expected!\n");
        exit(1);
    }
    if(fscanf(fp,"%d",&max_size)!=1){
        fprintf(stderr,"ERROR: maximum size expected!\n");
        exit(1);
    }
    if(fscanf(fp,"%d",&n_facs)!=1){
        fprintf(stderr,"ERROR: number of facilitites expected!\n");
        exit(1);
    }
    if(fscanf(fp,"%d",&n_clis)!=1){
        fprintf(stderr,"ERROR: number of clients expected!\n");
        exit(1);
    }

    // Initialize problem
    problem *prob = problem_init(n_facs,n_clis);
    prob->size_restriction_minimum = min_size;
    prob->size_restriction_maximum = max_size;
    prob->transport_cost = transport_cost;
    prob->client_gain = client_gain;
    prob->unassigned_cost = unassigned_cost;

    // Read facility costs
    for(int i=0;i<prob->n_facs;i++){
        if(fscanf(fp," %lf",&prob->facility_cost[i])!=1){
            fprintf(stderr,"ERROR: facility cost expected!\n");
            exit(1);
        }
    }
    // Read client weights and distances
    for(int i=0;i<prob->n_clis;i++){
        if(fscanf(fp,"%lf",&prob->client_weight[i])!=1){
            fprintf(stderr,"ERROR: client weight expected!\n");
            exit(1);
        }
        for(int j=0;j<prob->n_facs;j++){
            if(fscanf(fp," %lf",&prob->distance[j][i])!=1){
                fprintf(stderr,"ERROR: distance expected!\n");
                exit(1);
            }
        }
    }

    // Default parameters
    prob->filter = BETTER_THAN_SUBSETS;
    prob->lower_bound = -INFINITY;
    prob->branch_and_bound = 1;

    return prob;
}

problem *load_orlib_format(FILE *fp){

    // Read the number of facilities:
    int n_facs;
    if(fscanf(fp,"%d",&n_facs)!=1){
        fprintf(stderr,"ERROR: number of facilities expected!\n");
        exit(1);
    }

    // Read the number of clients:
    int n_clis;
    if(fscanf(fp,"%d",&n_clis)!=1){
        fprintf(stderr,"ERROR: number of clients expected!\n");
        exit(1);
    }

    problem *prob = problem_init(n_facs,n_clis);

    // For each facility
    int cap_0_warn = 0;
    for(int i=0;i<prob->n_facs;i++){

        // Read facility capacity
        double capacity = 0;
        char cap_text[200];
        if(fscanf(fp,"%s",cap_text)!=1){
            fprintf(stderr,"ERROR: facility capacity expected!\n");
            exit(1);
        }
        if(strcmp(cap_text,"capacity")!=0){
            if(sscanf(cap_text,"%lf",&capacity)!=1){
                fprintf(stderr,"ERROR: facility capacity isn't valid!\n");
                exit(1);
            }
        }
        if(capacity!=0 && !cap_0_warn){
            fprintf(stderr,"WARNING: facility capacity is %lf (not 0)! IGNORING.\n",capacity);
            cap_0_warn = 1;
        }

        // Read facility cost
        if(fscanf(fp,"%lf",&prob->facility_cost[i])!=1){
            fprintf(stderr,"ERROR: facility %d cost expected!\n",i);
            exit(1);
        }
    }

    // For each client
    int all_demands_0 = 1;
    for(int j=0;j<prob->n_clis;j++){
        // Read client demand
        double demand;
        if(fscanf(fp,"%lf",&demand)!=1){
            fprintf(stderr,"ERROR: client %d demand expected!\n",j);
        }
        if(demand!=0) all_demands_0 = 0;
        prob->client_weight[j] = demand;

        // Add distances to facility-city matrix:
        for(int i=0;i<prob->n_facs;i++){
            double dist;
            if(fscanf(fp,"%lf",&dist)!=1){
                fprintf(stderr,"ERROR: cost from facility %d to client %d expected!\n",i,j);
                exit(1);
            }
            if(demand==0){
                assert(all_demands_0 || dist==0);
                prob->distance[i][j] = dist;
            }else{
                prob->distance[i][j] = dist/demand;
            }
        }
    }

    // Unsetted values:
    prob->size_restriction_minimum = -1; // SPLP
    prob->size_restriction_maximum = -1; // SPLP
    prob->transport_cost = 1;
    if(all_demands_0){
        for(int j=0;j<prob->n_clis;j++){
            prob->client_weight[j] = 1;
        }
    }
    prob->client_gain = 0;
    prob->unassigned_cost = INFINITY;
    prob->filter = BETTER_THAN_SUBSETS;
    prob->lower_bound = -INFINITY;
    prob->branch_and_bound = 1;

    //
    return prob;
}

problem *new_problem_load(const char *file){
    printf("Reading file \"%s\"...\n",file);
    FILE *fp = fopen(file,"r");
    if(fp==NULL){
        fprintf(stderr,"ERROR: couldn't open file \"%s\"!\n",file);
        exit(1);
    }
    // Read first string to check if it is on SIMPLE format
    char buffer[400];
    if(fscanf(fp,"%s",buffer)!=1){
        fprintf(stderr,"ERROR: couldn't read first string!\n");
        exit(1);
    }

    fseek(fp,0,SEEK_SET); // Reset reading
    problem *prob;

    // Check if it is simple format
    if(strcmp(buffer,"FILE:")==0){
        printf("SIMPLE format identified.\n");
        prob = load_simple_format(fp);
    }else{
        buffer[2] = '\0';
        // Check if it is a DC format
        if(strcmp(buffer,"DC")==0){
            int dc_version = 0;
            int n_read = sscanf(&buffer[3],"V%d",&dc_version);
            if(n_read<1){
                fprintf(stderr,"ERROR: couldn't identify DC format version \"%s\"!\n",&buffer[3]);
                exit(1);
            }
            printf("DC format identified, version %d.\n",dc_version);
            prob = load_dc_format(fp,dc_version);
        }
        // Assume ORLIB format
        else{
            printf("ORLIB format assumed.\n");
            prob = load_orlib_format(fp);
        }
    }

    // Close file
    fclose(fp);
    printf("Done reading.\n");

    return prob;
}