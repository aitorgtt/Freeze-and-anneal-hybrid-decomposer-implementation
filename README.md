# Freeze and Anneal hybrid implementation.

Implementation of the iterative process of freezing variables in the Freeze and Anneal quantum annealing-based hybrid optimization solver. The idea behind this code was first introcuded in https://link.springer.com/article/10.1007/s42484-021-00039-9

It has been developed using D-Wave Ocean SDK python classes and Pygad for the genetic algorithm.

The implementation is in the form of a workflow similar to the ones offered by D-Wave in the dwave-hybrid python class, such as the EnergyImpactDecomposer or the QPUSubproblemAutoEmbeddingSampler.
