#include "load.h"
#include "redstrategy.h"
#include "problem.h"
#include "rundata.h"
#include "solution.h"

/*
The opt_checker is a simple program that reads the facilitites
present in the solution on a .opt or .bub file and uses those
indexes to build its own solution and retrieve it in the same
format.
*/

#include <time.h>
#include <sys/time.h>

int solution_check_not_optimal_assigns(const problem *prob, const solution *sol){
    int bad_assigns = 0;
    for (int i=0; i<prob->n_clis; i++){
        // Find the best value that client i could be assigned to
        int best_f = -1;
        double value_best = -INFINITY;
        for(int k=0;k<sol->n_facs;k++){
            int f = sol->facs[k];
            double value_alt = problem_assig_value(prob,f,i);
            if(best_f==-1 || value_alt>value_best){
                value_best = value_alt;
                best_f = f;
            }
        }
        // Current facility value:
        int c = sol->assigns[i];
        double value_c = problem_assig_value(prob,c,i);
        if(value_c < value_best){
            fprintf(stderr,"WARNING: Not optimal assign. %d->%d should be %d->%d.\n",
                c,i,best_f,i);
            bad_assigns += 1;
        }
    }
    return bad_assigns;
}

int main(int argc, const char **argv){
    // Print information if arguments are invalid
    if(argc!=3){
        fprintf(stderr,"usage: %s <input> <input_opt>\n",argv[0]);
        exit(1);
    }

    const char *input_fname = argv[1];
    const char *opt_fname = argv[2];

    // Read problem
    problem *prob = new_problem_load(input_fname);

    // Create empty solution
    solution *solution = solution_empty(prob);

    // Read opt file
    FILE *fp = fopen(opt_fname,"r");
    assert(fp!=NULL);
    for(int i=0;i<prob->n_clis;i++){
        int fac;
        int n_read = fscanf(fp,"%d",&fac);
        assert(n_read==1);
        solution_add(prob,solution,fac,NULL);
        // Get the assign cost by the assignment made by the algorithm
        double old_assign_cost = problem_assig_value(prob,solution->assigns[i],i);
        // Get the assign cost by the assignment on the .opt solution
        double opt_assign_cost = problem_assig_value(prob,fac,i);
        // Assert that both are the same
        assert(old_assign_cost==opt_assign_cost);
        // Replace algorithm assignment with .opt assignment
        solution->assigns[i] = fac;
    }
    fclose(fp);

    // Print result
    for(int i=0;i<prob->n_clis;i++){
        printf("%d ",solution->assigns[i]);
    }
    printf("%.9lf\n",-solution->value);

    // Perform one last additional check
    assert(solution_check_not_optimal_assigns(prob,solution)==0);

    // Free memory
    solution_free(solution);
    problem_free(prob);
}
