# BRKGA 

## Compile and Execute BRKGA

1. make
2. ./brkga-grundy --inst path_to_instance --solver solver_type --seed $seed --time 300 --pe 0.3 --pm 0.10 --rhoe 0.6 --factor 3 --reset 1 --ls 1 --nToImprove 5 --improvementType 2 --firstImprovement 1 --loop 1

### Possible values for solver_type

- grundyNP BRKGA decoder for a Grundy coloring.
- grundyCP BRKGA decoder for a connected Grundy coloring.
- grundyNPopt BRKGA decoder to be executed with local search for a Grundy coloring (run with --ls 1).
- grundyCPopt BRKGA decoder to be executed with local search for a connected Grundy coloring (run with --ls 1).

### Possible values for pe pm

- pe is the percentage of the elite population.
- pm is the percentage of the mutant population.
- Any value in the range (0,1)
- We indicate that "pe" is in the range 0.1 and 0.3.
- We indicate that "pm" is in the range 0.1 and 0.25.

### Possible values for rhoe

- Indicates the probability of inheritance.
- By definition "rhoe" > 0.5.

### Possible values for time

- Any integer value > 0. 
- Time is in seconds.
- We use the default timeout of 300 seconds.

### Possible values for factor

- Population size is defined as chromosome size times "factor".
- Any integer number > 0.
- The most common options are 1, 2, 3.

### Possible values for seed

- Any integer number

### Possible values to reset

- 1: use restart strategy after 2000 generations without found a new best solution.
- 0: do not use restart strategy.

### Possible values for ls

- 1: run local search
- 0: do not run local search

### Possible values to nToImprove

- Number of chromossomes to run the local search.
- Any integer number between 1 and population size.
- Recommendation: use the default value of 5 chromosomes.

### Possible values to improvementType

- 0: select the n best chromosomes.
- 1: select ceil(n/2) and floor(n/2) best chromosomes from elite and non-elite population, respectively.
- 2: select the best chromosome and (n-1) random chromosomes in the elite population.

### Possible values to firstImprovement

- 1: use first improvement strategy in the local search. 
- 0: use best improvement strategy in the local search.

### Possible values to loop

- 1: In local search, the search will continue from the selected improved solution.
- 0: In local search, it will terminate as soon as a new best solution is selected or it finishes exploring the neighborhood of the current solution without finding any better solution.