"""
This program generates an example problem in the DC_V1 format.
It locates clients and facilities in random positions in euclidean
then computes the distance between them and outputs a file.
"""

import sys
import numpy as np

def main():
    # Parse input arguments
    if len(sys.argv)!=3:
        print("usage: %s <n> <output>"%sys.argv[0])
        sys.exit(1);
    n = int(sys.argv[1])
    m = int(1+n*(np.random.random()*1.5+0.5))
    out_fname = sys.argv[2]
    # Client and facility positions
    dims = 3
    facs = np.random.random((n,dims))
    clis = np.random.random((m,dims))
    # Distances
    deltas = facs[np.newaxis,:,:] - clis[:,np.newaxis,:]
    dists = 100*np.sum(deltas**2,axis=-1)
    # Parameters
    C = 0.15
    P = 0.10
    client_weights = 2*np.random.random(m)
    client_gain = 100*(0.25+1.5*np.random.random(1))*(3*P/np.pi)**0.5
    transport_cost = 1*(0.25+1.5*np.random.random(1))
    unassigned_cost = 100*(0.1*np.random.random(1))
    facility_costs = 100*(0.25+1.5*np.random.random(n)*(C*client_gain*n*P))
    # Save the solution
    fi = open(out_fname,"w")
    fi.write("DC_V1\n")
    fi.write("%.3f\n"%transport_cost)
    fi.write("%.3f\n"%client_gain)
    fi.write("%.3f\n"%unassigned_cost)
    fi.write("%d %d\n"%(-1,-1))
    fi.write("%d %d\n"%(n,m))
    for i in range(n):
        fi.write(" %.2f"%facility_costs[i])
    fi.write("\n")
    for i in range(m):
        fi.write("%.2f\n"%client_weights[i])
        for j in range(n):
            fi.write(" %.2f"%dists[i][j])
        fi.write("\n")
    fi.close()

if __name__ == "__main__":
    main()