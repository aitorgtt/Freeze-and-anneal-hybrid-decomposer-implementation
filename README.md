# Freeze and Anneal hybrid implementation

Implementation of the iterative process of freezing variables in the Freeze and Anneal quantum annealing-based hybrid optimization solver. The idea behind this code was first introcuded in https://link.springer.com/article/10.1007/s42484-021-00039-9

QUBO problem decomposing method that takes the original problem and finds a subproblem of a given size (number of variables), which in theory, will be the most conflictive part of the original problem. The process to reduce the problem consists on running a genetic algorithm to solve the problem, to them freeze the variables which stayed constant on the 80% of all the solutions of the final population. Freezing variables means making them constant to the value they have taken on the mayority of solutions. This process is repeated until the problem gets reduced to the desired size.

It has been developed using D-Wave Ocean SDK python classes and Pygad for the genetic algorithm.

The implementation is in the form of a workflow similar to the ones offered by D-Wave in the dwave-hybrid python class, such as the EnergyImpactDecomposer or the QPUSubproblemAutoEmbeddingSampler.
