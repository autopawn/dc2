# dc2

**dc2** is a solver for uncapacitated facility location problems like SPLP and p-median, based on Disperse Construction, a constructive incomplete search method that explores the search space using diversity.

This is the second version of the original [dc_splp](https://github.com/autopawn/dc_splp).

It is mainly composed of:
* **Expansion**: Straightforward method to generate all possible solutions of size `n+1` by adding one facility to solutions size `n`.
* **Filtering**: A method that deletes some of the solutions if they are not worth exploring, i.e. leading to worse solutions that the one from where it originated.
* **Reduction**: A method to select a representative subset of solutions from a given set.
* **Local search**: A method that perform moves on a given solution until no further moves can increase the solution value.

![](imgs/dc2.png)

The search advances using these methods on several iterations, as shown in the previous diagram, until no more solutions are present. Additionaly it uses Branch and Bound to further filter solutions.

## Model

The solver aims to solve the following general uncapacitated location problem model, on its combinatorial form:

![](imgs/model_formula.gif)

where:

* `S` is a solution, the set of indexes of the opened facilities, a subset of `{0,1,2,...,n}`
* `C` is the set of all clients indexes `{0,1,2,...,m-1}`.
* `w_j` is the weight of client `j`.
* `K` is the cost of having clients unassigned.
* `C` is the profit for reaching clients.
* `T` is the transport cost for each unit of distance.
* `d_ij` is the distance between facility `i` and client `j`.
* `f_i` is the cost of opening the facility `i`.
* `s_min` is the minimum size that the solutions can have.
* `s_max` is the maximum size that the solutions can have.

SPLP, p-median, set covering, and maximum coverage problems can be **reduced** to this model:

| Problem | Values |
| :------ | ------ |
| SPLP | `K=inf`, `C=0` |
| p-median | `K=inf`, `C=0`, `f_i=0`, `s_max=p` |
| set covering | `K=n+1`, `C=0`, `f_i=1`, `T=inf`, `d_ij={0 or 1}` |
| maximum coverage | `K=0`, `C=1`, `f_i=0`, `T=inf`, `d_ij={0 or 1}`, `s_max=k` |

Furthermore, the solver can read a SPLP and a p-median problem file formats directly.

# Usage

Compile the program using `make`
```bash
make
```

Then execute it as follows:

```bash
./bin/dc [FLAGS] {reduction:k} <input> <output>
```

for example:
```
./bin/dc scripts/example_problem out.txt
```
solves the problem on `scripts/example_problem` and saves the solution to `out.txt`.

Normally the default parameters are suitable for getting good solutions, but they can be specified, for example:
```
./bin/dc -t8 -l rank:2000 sdce+:100 scripts/example_problem out.txt
```

executes the solver using 8 threads with a two step reduction method (first sampling based on `rank`, then `sdce+` until `100` solutions are reached) and skips local searches.

# Formats supported

The solver currently supports 3 formats, the ORLIB-cap format and the Simple format, specified on the [UflLib benchmark](https://resources.mpi-inf.mpg.de/departments/d1/projects/benchmarks/UflLib/data-format.html), and the more general DC_V1 format.

## ORLIB-cap format

Intended for SPLP problems.

The first line has the amount of facilities and clients of the problem:
```
[n] [m]
```
The next `n` lines consist on the opening cost and capacity each facility. NOTE: as this solver is meant for uncapacitated location problems, capacities are expected to be left at 0.
```
[capacity] [f_i]
```
The next `2m` lines, `2` for each client, contain:
* The demand (weight) of the client (if marked as 0, is set to 1).
* The cost of allocating all that demand (i.e. distance*weight) on each of the facilitites .
```
[w_j]
[d_0j*w_j] [d_1j*w_j] ... [d_(n-1)j*w_j]
```

Example: `n=4`, `m=3`:
```
4 3
0 300
0 400
0 150
0 200
0
130 140 130 100
0
120 100 90 120
0
80 50 140 150
```

## Simple format

Used for SPLP and p-median.

The first line consists on `FILE: ` and a filename.
```
FILE: [filename]
```
The next line contains the number of facilitites, clients, and size restriction:
```
[n] [m] [s_max]
```
The next `n` lines consist on the number of the facility (starting at 1), the opening costs, and the transport cost to each client
```
[i+1] [f_i] [d_i0] [d_i1] ... [d_i(m-1)]
```

Example: `n=4`, `m=3`:
```
FILE: Example.txt
4 3 0
1 300 130 120 80
2 400 140 100 50
3 150 130 90 140
4 200 100 120 150
```

## DC_V1 format

The DC_V1 format is more complex as it allows to set all model variables.

The fist line contains the string `DC_V1`, used to identify the format
```
DC_V1
```
The following lines contain the transport cost, the client gain and the unassigned cost multipliers.
```
[T]
[C]
[K]
```
The following line contains the minimum and maximum solutions sizes, can be left as `-1` if they don't apply.
```
[s_min] [s_max]
```
The following line contains the number of facilitites and clients.
```
[n] [m]
```
The following line contains the `n` facility costs
```
[f_0] [f_1] ... [f_(n-1)]
```
The following `2m`, `2` for each client contain:
* The client weight.
* The distances from this client to each facility.
```
[w_j]
[d_0j] [d_1j] ... [d_(n-1)j]
```

Example: `n=4`, `m=5`:
```
DC_V1
10.0
20.0
inf
-1 3
4 5
 200 400 300 400
1.0
 10 20 30 40
2.0
 15 5 5 15
3.0
 40 30 20 10
4.0
 25 15 5 5
3.0
 5 5 15 25
```

# Parameters

## Flags

The following flags can be used to specify different behaviours:

| Flag | Effect |
| :--- | ------ |
| `-b` | Don't perform Branch & Bound as additional filter. |
| `-l` | Skip local searches | 
| `-r<n>` | Sets the random seed to `n`, so execution is deterministic. |
| `-n<n>` | Sets the number of target solutions (1 by default). <br> Use with `-b` to get diverse solutions. |
| `-R<n>` | Performs `n` restarts, useful with random reduction components. <br> B&B bound is kept after restarts.  |
| `-t<n>` | The number of threads to use. |
| `-s<n>` | Sets the minimum size to `n`. <br> Solutions of smaller size are not considered as results. <br> Local search is not performed on them. |
| `-S<n>` | Sets the maximum size to `n`. <br> Once it is reached, the iteration stops.
| `-f<n>` | The filter level, can range from 0 to 4: <br> `-f0`: don't filter any solution. <br> `-f1`: solution should be better than the empty solution. <br> `-f2`: solution should be better than its worst parent. <br> `-f3`: solution should be better than its best parent. <br> `-f4`: solution should be better than any possible parent (default) | 

## Reduction strategies

The reduction is performed chaining one or more reduction strategies. These reduction strategies vary in the results they may have, memory and time complexity.

The default (and recommended) reduction strategy, is:
```
rand1:6000 sdce+:200
```
which picks 6000 solutions randomly and then applies scde to select 200.

**If you want to invest more computational power to solve the problem, you could indicate a reduction strategy that selects more solutions**:
```
rand1:10000 sdce+:400
```

You may also skip the previous random selection, which will be more costly but will result on more representative solutions.
```
sdce+:400
```

Complex reduction strategies like `sdce+` can work with a given dissimilitude metric, so for instance if your problem is metric, you may use `sdce+:400:mgemin` and `sdce+:400:mgesum` otherwise. These dissimilitudes require precomputations.

### Strategies

The **simple** strategies are very fast and do not use a dissimilitude metric between pairs of solutions, so are a good first step before using more complex strategies:

| **Strategy** | **Description** |
| :----------  | :-------------- |
| `rand:<n>`   | Pick `n` solutions at random. <br> Uniform probability. |
| `rand1:<n>`  | Same as `rand` but always pick the best solution.  |
| `rank:<n>`   | Pick `n` solutions at random. <br> Probability proportional to the reciprocal of the rank <br> (the position in a list sorted by value). |
| `rank1:<n>`  | Same as `rank` but always pick the best solution.  |
| `best:<n>`   | Pick the best `n` solutions. |

The **complex** strategies make use of a **dissimilitude** metric to compare between solutions and thus, pick an spacially different set of solutions.

| **Strategy** | **Description** |
| :----------  | --------------- |
| `sdce:<n>:<di>` | Select `n` solutions using Glover's simple diversity-based <br> clustering initialization method, enhanced. <br> Using the `<di>` dissimilitude metric (default: `pcd`).
| `sdce+:<n>:<di>` | Same as `scde` but the bests solutions of each cluster <br> are selected instead of the centroids.
| `vrh:<n>:<di>:<v>` | Select `n` solutions using the VR-Heuristic <br> with vision range `v` (default `2n`). <br> Using the `<di>` dissimilitude metric (default: `pcd`).

The following table lists the complexities to select `Q` solutions from a set of size `P`.

| **Strategy** | **Dissimilitudes** | **Memory** |
| :----------  | :-------------: | :--------: |
| `sdce` | `P Q` | `O(P)` |
| `sdce+` | `P Q` | `O(P)` |
| `vrh` | `2 P v` | `O(P v)` |

### Dissimilitude metrics:

The following dissimilitude metrics are available:

| `<di>` | **Description** |
| :---- | ------- |
| `mgemin` | Mean Geometric Error. <br> Using **min triangle** as facility-facility distance. |
| `mgesum` | Mean Geometric Error. <br> Using **sum of deltas** as facility-facility distance. |
| `haumin` | Hausdorff distance. <br> Using **min triangle** as facility-facility distance. | 
| `hausum` | Hausdorff distance. <br> Using **sum of deltas** as facility-facility distance. |
| `pcd`    | Per client delta. <br> Doesn't use facility-facility distances. |

Facility-facility distances:
* **min triangle**:
    ```
    df(a,b) = min_j d(a,j)+d(b,j)
    ```
    This facility-facility distance is the best for **metric** problems, but may be very bad for **non-metric** problems. Its intended to be geographical distance.
    
    Precomputation of all them costs O(n^2 m).

* **sum of deltas**:
    ```
    df(a,b) = sum_j |v(a,j)-v(b,j)|
    ```
    This facility-facility distance is the best for **non-metric** problems. Compares the assignment costs.

    Precomputation of all them costs O(n^2 m).

Solution-solution dissimilitudes:
* **Mean geometric error**:
    ```
    D(A,B) = sum_a inf_b df(a,b) + sum_b inf_a df(b,a)
    ```
    Is stable and brings good results, however it costs is proportional to O(p^2) where p is the size of the solutions.
* **Hausdorff**:
    ```
    D(A,B) = max {sup_a inf_b df(a,b), sup_b inf_a df(b,a)}
    ```
    Is faster but may have poor results, it is prone to ties.
* **Per client delta**:
    ```
    D(A,B) = sum_j |v(A,j)-v(B,j)|
    ```
    Uses the client assignment costs as a vector which is compared, its cost is O(m) (number of clients) so may be good for problems with very large solutions.

## Local search

The local searches are performed using Whitaker's fast swap heuristic.
