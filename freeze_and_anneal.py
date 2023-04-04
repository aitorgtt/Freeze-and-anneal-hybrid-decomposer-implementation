import hybrid
from hybrid.utils import bqm_induced_by, updated_sample
from hybrid.core import State
import pygad
import numpy as np


class FreezeAnnealDecomposer(hybrid.Runnable):

    def __init__(self, size, num_generations, num_parents_mating, proportion=0.2, fitness_batch_size=None, initial_population=None,
                 sol_per_pop=None, num_genes=None, init_range_low=-4, init_range_high=4, parent_selection_type="sss",
                 keep_parents=-1, keep_elitism=1, K_tournament=3, crossover_type="single_point", crossover_probability=None,
                 mutation_type="random",mutation_probability=None, mutation_by_replacement=False, mutation_percent_genes='default',
                 mutation_num_genes=None,  random_mutation_min_val=-1.0, random_mutation_max_val=1.0,
                 allow_duplicate_genes=True, on_start=None, on_fitness=None, on_parents=None, on_crossover=None,
                 on_mutation=None, callback_generation=None, on_generation=None, on_stop=None, delay_after_gen=0.0,
                 save_best_solutions=False, save_solutions=False,
                 suppress_warnings=False,
                 stop_criteria=None,
                 parallel_processing=None,
                 random_seed=None, **runopts):
        self.size = size
        self.proportion = proportion
        self.num_generations = num_generations
        self.num_parents_mating = num_parents_mating
        self.fitness_batch_size = fitness_batch_size
        self.initial_population = initial_population
        self.sol_per_pop = sol_per_pop
        self.num_genes = num_genes
        self.init_range_low = init_range_low
        self.init_range_high = init_range_high
        self.parent_selection_type = parent_selection_type
        self.keep_parents = keep_parents
        self.keep_elitism = keep_elitism
        self.K_tournament = K_tournament
        self.crossover_type = crossover_type
        self.crossover_probability = crossover_probability
        self.mutation_type = mutation_type
        self.mutation_probability = mutation_probability
        self.mutation_by_replacement = mutation_by_replacement
        self.mutation_percent_genes = mutation_percent_genes
        self.mutation_num_genes = mutation_num_genes
        self.random_mutation_min_val = random_mutation_min_val
        self.random_mutation_max_val = random_mutation_max_val
        self.allow_duplicate_genes = allow_duplicate_genes
        self.on_start = on_start
        self.on_fitness = on_fitness
        self.on_parents = on_parents
        self.on_crossover = on_crossover
        self.on_mutation = on_mutation
        self.callback_generation = callback_generation
        self.on_generation = on_generation
        self.on_stop = on_stop
        self.delay_after_gen = delay_after_gen
        self.save_best_solutions = save_best_solutions
        self.save_solutions = save_solutions
        self.suppress_warnings = suppress_warnings
        self.stop_criteria = stop_criteria
        self.parallel_processing = parallel_processing
        self.random_seed = random_seed
        super(FreezeAnnealDecomposer, self).__init__(**runopts)

    @staticmethod
    def _bqm_to_tabu_qubo(bqm):
        # construct dense matrix representation
        ldata, (irow, icol, qdata), offset, varorder = bqm.binary.to_numpy_vectors(return_labels=True)
        ud = np.zeros((len(bqm), len(bqm)), dtype=np.double)
        ud[np.diag_indices(len(bqm), 2)] = ldata
        ud[irow, icol] = qdata
        return ud


    def fitness_func(self, solution, solution_idx):
        fitness = - np.dot(solution,np.matmul(self.qubo, solution))
        return fitness

    def min_variance_index(self, num_genes, sol_per_pop, newsampleset, best_solution):
        # Search for the variable with the highest variance in the last population
        prop = np.zeros(num_genes)
        for var in range(num_genes):
            for sol in range(sol_per_pop):
                prop[var] += np.abs(newsampleset[sol][var] - best_solution[var])
            prop[var] = prop[var]/sol_per_pop
        best = np.argmin(prop)
        if prop[best] < self.proportion:
            return best
        return None

    def next(self, state, **runopts):
        state = state.updated(subproblem=state.problem) # Create subproblem attribute
        n = len(state.samples.first.sample) - self.size
        variables = list(range(len(state.samples.first.sample)))

        for ite in range(n):
            bqm = state.subproblem
            self.qubo = self._bqm_to_tabu_qubo(bqm)
            fitness_function = self.fitness_func
            num_genes = len(bqm)
            ga_instance = pygad.GA(num_generations = self.num_generations, num_parents_mating = self.num_parents_mating,
                 fitness_func = fitness_function, num_genes = num_genes, save_solutions=True, fitness_batch_size=self.fitness_batch_size,
                 initial_population=self.initial_population, sol_per_pop=self.sol_per_pop, init_range_low=self.init_range_low,
                 init_range_high=self.init_range_high, gene_type=int, parent_selection_type=self.parent_selection_type, keep_parents=self.keep_parents,
                 keep_elitism=self.keep_elitism, K_tournament=self.K_tournament, crossover_type=self.crossover_type, crossover_probability=self.crossover_probability,
                 mutation_type=self.mutation_type, mutation_probability=self.mutation_probability, mutation_by_replacement=self.mutation_by_replacement,
                 mutation_percent_genes=self.mutation_percent_genes, mutation_num_genes=self.mutation_num_genes, random_mutation_min_val=self.random_mutation_min_val,
                 random_mutation_max_val=self.random_mutation_max_val, gene_space=[0,1], allow_duplicate_genes=self.allow_duplicate_genes,
                 on_start=self.on_start, on_fitness=self.on_fitness, on_parents=self.on_parents, on_crossover=self.on_crossover, on_mutation=self.on_mutation,
                 callback_generation=self.callback_generation, on_generation=self.on_generation, on_stop=self.on_stop, delay_after_gen=self.delay_after_gen,
                 save_best_solutions=self.save_best_solutions, suppress_warnings=self.suppress_warnings, stop_criteria=self.stop_criteria,
                 parallel_processing=self.parallel_processing, random_seed=self.random_seed)
            ga_instance.run()
            best_solution = ga_instance.best_solution()[0]

            # Save the latest solution in the sample attribute
            sample = state.samples.change_vartype(bqm.vartype).first.sample
            new_sample = {variables[i]: best_solution[i] for i in range(num_genes)}
            state = State.from_sample(updated_sample(sample, new_sample), state.problem)

            print(state.samples.first.energy)

            # Create subproblem freezing the variable with the lowest variance, is this stays constant al least the self.proportion% of the times
            newsampleset = ga_instance.solutions[-self.sol_per_pop:]
            freeze_index = self.min_variance_index(num_genes, self.sol_per_pop, newsampleset, best_solution)
            if freeze_index is None:
                print("All variables change too much, no possible freezing. Try with another method.")
                break
            variables.pop(freeze_index)
            subbqm = bqm_induced_by(bqm, set(variables), sample)
            state = state.updated(subproblem=subbqm)

        return state
