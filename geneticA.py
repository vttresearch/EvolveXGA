import random
import copy
import numpy as np
from operator import attrgetter
import multiprocessing as mp
from joblib import Parallel, delayed

from six.moves import range
import itertools
from itertools import combinations


class GeneticAlgorithm(object):
    """Genetic Algorithm class.
    This is the main class that controls the functionality of the Genetic
    Algorithm.
    A simple example of usage:
    >>> # Select only two items from the list and maximise profit
    >>> from pyeasyga.pyeasyga import GeneticAlgorithm
    >>> input_data = [('pear', 50), ('apple', 35), ('banana', 40)]
    >>> easyga = GeneticAlgorithm(input_data)
    >>> def fitness (member, data):
    >>>     return sum([profit for (selected, (fruit, profit)) in
    >>>                 zip(member, data) if selected and
    >>>                 member.count(1) == 2])
    >>> easyga.fitness_function = fitness
    >>> easyga.run()
    >>> print easyga.best_individual()
    """

    def __init__(self,
                 seed_data, #this is now a single vector
                 testCompStart,
                 aminoAcids,
                 population_size=200,
                 generations=100,
                 crossover_probability=0.8,
                 mutation_probability=0.8, 
                 elitism=True,
                 maximise_fitness=True):
        """Instantiate the Genetic Algorithm.
        :param seed_data: input data to the Genetic Algorithm
        :type seed_data: list of objects
        :param int population_size: size of population
        :param int generations: number of generations to evolve
        :param float crossover_probability: probability of crossover operation
        :param float mutation_probability: probability of mutation operation
        """

        self.seed_data = seed_data
        self.testCompStart = testCompStart
        self.aminoAcids = aminoAcids
        self.population_size = population_size
        self.generations = generations
        self.crossover_probability = crossover_probability
        self.mutation_probability = mutation_probability
        self.elitism = elitism
        self.maximise_fitness = maximise_fitness

        self.current_generation = []
        self.scoreslist = []

        def create_individual(seed_data):
            """Create a candidate solution representation.
            e.g. for a bit array representation:
            >>> return [random.randint(0, 1) for _ in range(len(data))]
            :param seed_data: input data to the Genetic Algorithm
            :type seed_data: list of objects
            :returns: candidate solution representation as a list
            """
            chromosome = np.copy(seed_data)


            checking = 0
            while checking == 0:
                compoundAmount = 0
                compounds = []
                while compoundAmount < 3:
                    idx = random.randint(0, len(testCompStart)-1)  
                    if chromosome[idx] == 1:
                        compoundAmount = compoundAmount + 1
                        chromosome[idx] = 0 
                        compounds.append(testCompStart[idx])
                aminoAcidAmount = 0
                for i in compounds:
                    if i in aminoAcids:
                        aminoAcidAmount = aminoAcidAmount + 1
                if aminoAcidAmount == 3:
                    chromosome = np.copy(seed_data)
                    checking = 0
                else:
                    checking = 1


            idxAmount = 0
            while idxAmount < 4:
                idx = random.randint(len(testCompStart), len(seed_data) - 1)  
                if chromosome[idx] == 1:
                    idxAmount = idxAmount + 1
                    chromosome[idx] = 0    

            return chromosome
            

        def crossover(parent_1, parent_2):
            """Crossover (mate) two parents to produce two children.
            :param parent_1: candidate solution representation (list)
            :param parent_2: candidate solution representation (list)
            :returns: tuple containing two children
            """
            
            indexes_parent_1 = np.argwhere(parent_1.genes==0).tolist()
            indexes_parent_2 = np.argwhere(parent_2.genes==0).tolist()
            reactionsToCrossover = []
            compoundsToCrossover = []
            for rxn in indexes_parent_1[0:3]:
                compoundsToCrossover.append(rxn)
            for rxn in indexes_parent_2[0:3]:
                if rxn not in compoundsToCrossover:
                    compoundsToCrossover.append(rxn)
    
            for rxn in indexes_parent_1[3:7]:
                reactionsToCrossover.append(rxn)
            for rxn in indexes_parent_2[3:7]:
                if rxn not in reactionsToCrossover:
                    reactionsToCrossover.append(rxn)
                
            
            combCompounds = [x for x in itertools.combinations(compoundsToCrossover, 3)]
            checking = 0
            while checking == 0:
                compounds = []
                set_1_c = random.choice(combCompounds)
                set_1_c_list2d = list(set_1_c)
                set_1_c_list = list(itertools.chain.from_iterable(set_1_c_list2d))
                for i in set_1_c_list:
                    compounds.append(testCompStart[i])
                aminoAcidAmount = 0
                for i in compounds:
                    if i in aminoAcids:
                        aminoAcidAmount = aminoAcidAmount + 1
                if aminoAcidAmount == 3:
                    combCompounds.remove(set_1_c)
                    checking = 0
                else:
                    checking = 1
            combCompounds.remove(set_1_c)
            try:
                checking = 0
                while checking == 0:
                    compounds = []
                    set_2_c = random.choice(combCompounds)
                    set_2_c_list2d = list(set_2_c)
                    set_2_c_list = list(itertools.chain.from_iterable(set_2_c_list2d))
                    for i in set_2_c_list:
                        compounds.append(testCompStart[i])
                    aminoAcidAmount = 0
                    for i in compounds:
                        if i in aminoAcids:
                            aminoAcidAmount = aminoAcidAmount + 1
                    if aminoAcidAmount == 3:
                        combCompounds.remove(set_2_c)
                        checking = 0
                    else:
                        checking = 1
            except IndexError:  
                set_2_c = set_1_c
            
            
            combReactions = [x for x in itertools.combinations(reactionsToCrossover, 4)]
            set_1_r = random.choice(combReactions)
            combReactions.remove(set_1_r)
            try:
                set_2_r = random.choice(combReactions)
            except IndexError:  
                set_2_r = set_1_r

            for i in indexes_parent_1:
                parent_1.genes[i] = 1

            for i in set_1_c:
                parent_1.genes[i] = 0

            for i in set_1_r:
                parent_1.genes[i] = 0

            for i in indexes_parent_2:
                parent_2.genes[i] = 1

            for i in set_2_c:
                parent_2.genes[i] = 0

            for i in set_2_r:
                parent_2.genes[i] = 0

            return parent_1, parent_2
            

        def mutate(individual):
            """Reverse the bit of a random index in an individual."""
            indexes = np.argwhere(individual.genes==0)                
            mutate_index_c = random.randrange(0,3)
            mutate_index_r = random.randrange(3,7)
            compoundToMutate = indexes[mutate_index_c]
            reactionToMutate = indexes[mutate_index_r]
            individual.genes[compoundToMutate] = 1
            individual.genes[reactionToMutate] = 1

            control = 0
            while control == 0:
                mutateAmountCompounds = 0
                while mutateAmountCompounds < 1:
                    mutate_index = random.randrange(len(testCompStart))
                    if individual.genes[mutate_index] == 1:
                        mutateAmountCompounds = mutateAmountCompounds + 1
                        individual.genes[mutate_index] = 0
                indexes2d = np.argwhere(individual.genes==0).tolist()
                indexes = list(itertools.chain.from_iterable(indexes2d))
 
                compounds = []
                for i in indexes[0:3]:
                    compounds.append(testCompStart[i])
                aminoAcidAmount = 0
                for i in compounds:
                    if i in aminoAcids:
                        aminoAcidAmount = aminoAcidAmount + 1
                if aminoAcidAmount == 3:
                    individual.genes[mutate_index] = 1
                    control = 0
                else:
                    control = 1
   
            mutateAmountReactions = 0
            while mutateAmountReactions < 1:
                mutate_index =len(testCompStart) + random.randrange(len(individual.genes)-len(testCompStart))
                if individual.genes[mutate_index] == 1:    
                    mutateAmountReactions = mutateAmountReactions + 1
                    individual.genes[mutate_index] = 0   


        def random_selection(population):
            """Select and return a random member of the population."""
            return random.choice(population)

        def tournament_selection(population):
            """Select a random number of individuals from the population and
            return the fittest member of them all.
            """
            if self.tournament_size == 0:
                self.tournament_size = 2
            members = random.sample(population, self.tournament_size)
            members.sort(
                key=attrgetter('fitness'), reverse=self.maximise_fitness)
            return members[0]

        self.fitness_function = None
        self.tournament_selection = tournament_selection
        self.tournament_size = self.population_size // 10
        self.random_selection = random_selection
        self.create_individual = create_individual
        self.crossover_function = crossover
        self.mutate_function = mutate
        self.selection_function = self.tournament_selection

    def create_initial_population(self):
        """Create members of the first population randomly.
        """
        initial_population = []
        for _ in range(self.population_size):
            genes = self.create_individual(self.seed_data)
            individual = Chromosome(genes) #an individual is an instance of class Chromosome 
            initial_population.append(individual)
        self.current_generation = initial_population

    def calculate_population_fitness(self):
        """Calculate the fitness of every member of the given population using
        the supplied fitness_function.
        """
        results = {}
        queue = mp.Queue()
        queue.put(results)

        id = 0
        children = []
        for individual in self.current_generation:
            individual.id = str(id)
            p = mp.Process(target=self.fitness_function, args=(individual, self.seed_data, queue)) 
            children.append(p)
            id = id + 1
            p.start()
        
        for p in reversed(children):
            p.join()

        results = queue.get()
        scorecounting =[]
        for i in range(id):
            score = results[str(i)]
            scorecounting.append(score)
            for individual in self.current_generation:
                if individual.id == str(i):
                    individual.fitness =score
        self.scoreslist.append(max(scorecounting))


    def rank_population(self):
        """Sort the population by fitness according to the order defined by
        maximise_fitness.
        """
        self.current_generation.sort(
            key=attrgetter('fitness'), reverse=self.maximise_fitness)


    def create_new_population(self):
        """Create a new population using the genetic operators (selection,
        crossover, and mutation) supplied.
        """
        new_population = []
        elite = copy.deepcopy(self.current_generation[0:10])
        selection = self.selection_function

        if self.elitism:
            new_population = elite
        

        while len(new_population) < self.population_size:
            parent_1 = copy.deepcopy(selection(self.current_generation))
            parent_2 = copy.deepcopy(selection(self.current_generation))

            child_1, child_2 = parent_1, parent_2
            child_1.fitness, child_2.fitness = 0, 0

            can_crossover = random.random() < self.crossover_probability
            can_mutate = random.random() < self.mutation_probability

            if can_crossover:
                child_1, child_2 = self.crossover_function(
                    parent_1, parent_2)
                

            if can_mutate:
                self.mutate_function(child_1)
                self.mutate_function(child_2)

        
            new_population.append(child_1)


            if len(new_population) < self.population_size:  
                new_population.append(child_2)



        self.current_generation = new_population

    def change_population(self):
        
        new_population = []
        elite = copy.deepcopy(self.current_generation[0])

        if self.elitism:
            new_population = [elite]

        while len(new_population) < self.population_size:
            genes = self.create_individual(self.seed_data)
            individual = Chromosome(genes) #an individual is an instance of class Chromosome 
            new_population.append(individual)

        self.current_generation = new_population
            
           

    def create_first_generation(self):
        """Create the first population, calculate the population's fitness and
        rank the population by fitness according to the order specified.
        """
        self.create_initial_population()
        self.calculate_population_fitness()
        self.rank_population()

    def create_next_generation(self):
        """Create subsequent populations, calculate the population fitness and
        rank the population by fitness in the order specified.
        """
        self.create_new_population()
        self.calculate_population_fitness()
        self.rank_population()

    def create_new_generation(self):
       """Creating new population every 100 generations so that it includes the 10 elite individuals from the previous generation, but the rest of the chromosomes are totally new.
       and ranking the poopultaion by fitness"""

       self.change_population()
       self.calculate_population_fitness()
       self.rank_population()


    def run(self,testCompStart,names):
        """Run (solve) the Genetic Algorithm."""
        file = open("best_individuals.csv","w+")
        self.create_first_generation()
        generationNumber = 1
        fitness, genes = self.best_individual()
        file.write(str(generationNumber) +"\n")
        file.write(str(fitness) +"\n")
        for i in range (len(genes)):
            if genes[i] == 0:
                if i < len(testCompStart):    
                    file.write(str(i)+","+ str(testCompStart[i]) +"\n")
                else:
                    file.write(str(i)+","+ str(names[i-len(testCompStart)]) +"\n")        
        for _ in range(1, self.generations):
            if _ % 2000 == 0:
                self.create_new_generation()
                generationNumber = generationNumber + 1
                fitness, genes = self.best_individual()
                file.write(str(generationNumber) +"\n")
                file.write(str(fitness) +"\n")
                for i in range (len(genes)):
                    if genes[i] == 0:
                        if i < len(testCompStart):    
                            file.write(str(i)+","+ str(testCompStart[i]) +"\n")
                        else:
                            file.write(str(i)+","+ str(names[i-len(testCompStart)]) +"\n")

            else:
                self.create_next_generation()
                generationNumber = generationNumber + 1
                fitness, genes = self.best_individual()
                file.write(str(generationNumber) +"\n")
                file.write(str(fitness) +"\n")
                for i in range (len(genes)):
                    if genes[i] == 0:
                        if i < len(testCompStart):    
                            file.write(str(i)+","+ str(testCompStart[i]) +"\n")
                        else:
                            file.write(str(i)+","+ str(names[i-len(testCompStart)]) +"\n")
        file.close()  

    def best_individual(self):
        """Return the individual with the best fitness in the current
        generation.
        """
        best = self.current_generation[0]
        return (best.fitness, best.genes)
    
    def best_individualGenes(self):
        best = self.current_generation[0]
        return (best.genes)

    def last_generation(self):
        """Return members of the last generation as a generator function."""
        return ((member.fitness, member.genes) for member
                in self.current_generation)


class Chromosome(object):
    """ Chromosome class that encapsulates an individual's fitness and solution
    representation.
    """

    def __init__(self, genes):
        """Initialise the Chromosome."""
        self.id = ""
        self.genes = genes
        self.fitness = -10000

    def __repr__(self):
        """Return initialised Chromosome representation in human readable form.
        """
        return repr((self.fitness, self.genes))







