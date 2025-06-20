import numpy as np
import cobra
from cobra.util import solver
from cobra import Model, Reaction, Metabolite
from math import floor
import cplex
import pandas as pd
import random
import geneticA
import matplotlib.pyplot as plt

def setBounds(cbmodel):
    for rxn in cbmodel.reactions:
        rxn.lower_bound = rxn.lower_bound * 10
        rxn.upper_bound = rxn.upper_bound * 10


def createSecretion(model, rxn):
    irxn = model.reactions.get_by_id(rxn)
    reaction = Reaction("secretion_" + rxn)
    reaction.name="secretion of rxn"
    reaction.lower_bound = 0
    reaction.upper_bound = 10000
    reaction.add_metabolites(irxn.metabolites)
    model.add_reactions([reaction])


def evolveX(comp, inh, model, targets, targetDirs, aer):
    #input:
    #comp -> components in the chemical environment, uptake rxn indeces
    #inh -> inhibitors in the chemical environment, indeces in the inhibition matrix
    #model -> Yeast 7.6 metabolic model in cobra toolbox format
    #targets -> vector of target flux indeces in the model
    #targetDirs -> a vector of target flux desired directions of change (UP/DOWN)
    #output:
    #strength -> the worst case of selection pressure on the target fluxes created by the individual chemical environment
    #coverage -> the optimal coverage obtainable under the above worst case selection pressure

    #mapping the chemical environment into the bounds of the model
    if aer == 1:
        #in practise no constraint on oxygen uptake or respiration
        model.reactions.get_by_id('r_0439').upper_bound =10000
    else:
        #no respiration, but oxygen for biosynthetic needs
        model.reactions.get_by_id('r_0439').upper_bound= 0


    #setting the bounds of glucose and ammonium uptake to 0
    model.reactions.get_by_id('r_1714').lower_bound = 0
    model.reactions.get_by_id('r_1654').lower_bound = 0


    #setting to 0 reactions from the model that should not be there
    model.reactions.get_by_id('r_4485').lower_bound = 0
    model.reactions.get_by_id('r_4485').upper_bound = 0
    model.reactions.get_by_id('r_4486').lower_bound = 0
    model.reactions.get_by_id('r_4486').upper_bound = 0
    model.reactions.get_by_id('r_4487').lower_bound = 0
    model.reactions.get_by_id('r_4487').upper_bound = 0
    model.reactions.get_by_id('r_4488').lower_bound = 0
    model.reactions.get_by_id('r_4488').upper_bound = 0


    #compounds in the chemical environment, allow uptake
    for rxn in comp:
        model.reactions.get_by_id(rxn).lower_bound = -10000
        model.reactions.get_by_id(rxn).upper_bound = 0
        #create the separate secretion reactions for the nutrients available for uptake
        createSecretion(model, rxn)


    #inhibitors in the chemical environment

    if len(inh) != 0:
        # find out whether they are annotated to genes or rxns
        genes = getGenes(model)
        if inh[1] in genes:
            # find the orfs affected, then identify rxns affected
            g = [x for x in genes if x not in inh]
            active_reactions = model.evaluate_gprs(g)
            for reaction in model.reactions:
                if reaction not in active_reactions:
                    model.reactions[reaction].set_lower_bound(0)
                    model.reactions[reaction].set_upper_bound(0)
        elif inh[1] in model2.reactions.keys():
            for reaction in model.reactions:
                if reaction not in inh:
                    model.reactions[reaction].set_lower_bound(0)
                    model.reactions[reaction].set_upper_bound(0)
        else:
            print('Incorrect inhibitor annotation matrix! Inhibitors not considered')

    #growth set to 1
    model.reactions.get_by_id('r_2111').lower_bound =1
    model.reactions.get_by_id('r_2111').upper_bound =1

    obj = {}
    for rxn in model.reactions:
        if rxn.id in comp:
            obj[rxn.reverse_variable] = -1
            #print(rxn)
    obj[model.reactions.get_by_id('r_2111').reverse_variable] = 0
    obj[model.reactions.get_by_id('r_2111').forward_variable] = 0	

    model.objective.set_linear_coefficients(obj)
    solution = model.optimize()

    #restricting nutrient uptake
    fluxsum= sum(solution.fluxes[comp])
    

    if len(solution.fluxes) > 0 and fluxsum > -75:
        tol = 1e-9
        problem = cplex.Cplex()
        problem.parameters.simplex.tolerances.optimality = tol
        # Creating the Aeq matrix

        # Adding variables to the problem (each reaction is a variable)
        for rxn in model.reactions:
            # the variables in targets with lower bound = 0 if they are UP get 1 in the objective and if they are DOWN -1
            # the rest of variables (targets negative lower bound and non targets) get 0 in the objective.
            if rxn.id in targets and rxn.lower_bound >= 0:
                ind = targets.index(rxn.id)
                dir = targetDirs[ind]
                if 'UP' in dir:
                    problem.variables.add(obj=[1], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])
                else:
                    problem.variables.add(obj=[-1], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])
            else:
                problem.variables.add(obj=[0], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])

        # Adding constrains to the problem (the stoichimetries of each metabolite in each reaction are the constrains)
        A = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame")
        names = A.index
        A = A.transpose()
        for name in names:
            r = A[A[name] != 0]
            problem.linear_constraints.add(lin_expr=[[r.index, r[name]]],
                                           senses=['E'],
                                           rhs=[0])



        #fluxes as variables
        for i in range(len(targets)):
            rxnTarget = model.reactions.get_by_id(targets[i])
            if 'UP' in targetDirs[i]:
                # if up-regulation target
                # the aim is to minimize the absolute values of the UP-regulation target
                # fluxes given optimal conversion of nutrients to biomass (=worst case scenario)
                if model.reactions.get_by_id(targets[i]).lower_bound <0:
                    # new variable
                    ubs = max(abs(rxnTarget.upper_bound),abs(rxnTarget.lower_bound))

                    problem.variables.add(obj=[1], lb=[0], ub=[ubs], names=["up_"+targets[i]])

                    const1 = [[targets[i], "up_"+targets[i]],[1,-1]]
                    const2 = [[targets[i], "up_"+targets[i]],[-1,-1]]
                    problem.linear_constraints.add(lin_expr=[const1, const2],
                                                       senses=['L','L'],
                                                       rhs=[0,0])



            elif 'DOWN' in targetDirs[i]:
                # if down - regulation target
                # the aim is to maximize the absolute values of the DOWN -
                # regulation target fluxes given optimal conversion of nutrients to biomass (=worst case scenarion)

                if model.reactions.get_by_id(targets[i]).lower_bound <0:
                    # 2 new variables
                    ubs = max(abs(rxnTarget.upper_bound),abs(rxnTarget.lower_bound))
                    problem.variables.add(obj=[-1], lb=[0], ub=[ubs], names=["down_" + targets[i]])
                    problem.variables.add(obj=[0], lb=[0], ub=[1], names=["down2_" + targets[i]])
                    problem.variables.set_types("down2_" + targets[i], problem.variables.type.binary)

                    const1 = [[targets[i], "down_" + targets[i]],[1, -1]]
                    const2 = [[targets[i], "down_" + targets[i]],[-1, -1]]
                    const3 = [[targets[i], "down_" + targets[i], "down2_" + targets[i]],[-1, 1, -20000]]
                    const4 = [[targets[i], "down_" + targets[i],"down2_" + targets[i]],[1, 1, 20000]]
                    problem.linear_constraints.add(lin_expr=[const1, const2, const3, const4],
                                                       senses=['L', 'L', 'L', 'L'],
                                                       rhs=[0, 0,0,20000])




        ub=-floor(solution.objective_value/tol)*tol
        cn = []
        coeff = []
        for rxn in comp:
            coeff.append(-1)
            cn.append(rxn)

        const1 = [cn, coeff]
        problem.linear_constraints.add(lin_expr=[const1],
                                       senses=['L'],
                                       rhs=[ub])

        problem.objective.set_sense(problem.objective.sense.minimize)
        
        try:
            problem.solve()
            value = -1000000
            if problem.solution.get_status() == 1:
                value = problem.solution.get_objective_value()
            return value

        except:
            return -1000000

    else:
        return -1000000





print("loading the model")

modelInit =cobra.io.read_sbml_model("Oxalate_pathway.xml")
modelInit.solver ='cplex'

model = modelInit.copy()
model.objective = "r_5020"

print("setting up")
testComp = []
df = pd.read_csv('uptake_ids_reduced.csv', squeeze=True) #here list of possible nutrients in evolution environment
testCompStart = df.values.tolist()
df2 = pd.read_csv('amino_acids.csv', squeeze=True)  
aminoAcids = df2.values.tolist()
inh =[]


targets =["r_5000", "r_5001", "r_5002"]
dir =["UP", "UP", "UP"]

print("evolveX")






# define reactions that we don't want the algorithm to suggest deletions to. These include essential reactions both for production and growth, reactions that take part in lipid metabolism, reactions in which there is no flux, exchange and transport reactions and reactions that don't have gene annotation..
from cobra.flux_analysis import flux_variability_analysis

#essential reactions for heterologous pathway

flux = flux_variability_analysis(model, fraction_of_optimum = 0.7)
pos=flux.loc[((flux['minimum']>0.000001) & (flux['maximum']>0.000001))]
neg=flux.loc[((flux['minimum']<-0.000001) & (flux['maximum']<-0.000001))]
result = pd.concat([pos,neg])
essentialProduction = result.index.tolist()

#essential reactions for growth

flux = flux_variability_analysis(modelInit, fraction_of_optimum = 0.1)
pos=flux.loc[((flux['minimum']>0.000001) & (flux['maximum']>0.000001))]
neg=flux.loc[((flux['minimum']<-0.000001) & (flux['maximum']<-0.000001))]
result = pd.concat([pos,neg])
essentialGrowth = result.index.tolist()

# reactions were there is no flux

flux = flux_variability_analysis(modelInit, fraction_of_optimum = 0.1)
zero=flux.loc[((flux['minimum']==0) & (flux['maximum']==0))]
result = pd.concat([zero])
zero = result.index.tolist()

#reactions that take part in lipid metabolism

lipiddf = pd.read_csv('lipids.csv', squeeze=True)
lipids = lipiddf.values.tolist()

#transport reactions

transport =[]
for rxn in modelInit.reactions:
    if 'transport' in rxn.name:
        transport.append(rxn.id)

#exchange reactions

exchange =[]
for rxn in modelInit.reactions:
    if 'exchange' in rxn.name:
        exchange.append(rxn.id)

#diffusion reactions

diffusion = []
for rxn in modelInit.reactions:
    if 'diffusion' in rxn.name:
        diffusion.append(rxn.id)

#reactions without gene annotation

noGeneAnnotation = []
for rxn in modelInit.reactions:
    gene = rxn.gene_reaction_rule
    if 'Y' not in gene:
        noGeneAnnotation.append(rxn.id)

#more reactions to remove

moreReactionsToRemovedf = pd.read_csv('more_reactions_to_remove.csv', squeeze=True)  #here list of reactions that should be removed in addition
moreReactionsToRemove = moreReactionsToRemovedf.values.tolist()

reactions = []
for rxn in modelInit.reactions:
    reactions.append(rxn.id)

for rxn in essentialProduction:
        reactions.remove(rxn)

for rxn in essentialGrowth:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in zero:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in lipids:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in transport:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in exchange:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in diffusion:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in noGeneAnnotation:
    if rxn in reactions:
        reactions.remove(rxn)

for rxn in moreReactionsToRemove:
    if rxn in reactions:
        reactions.remove(rxn)


reactionsToBeRemoved = []
for rxn in modelInit.reactions:
    if rxn.id not in reactions:
        reactionsToBeRemoved.append(rxn)

data = []
data = np.ones(len(testCompStart)+len(modelInit.reactions)-len(reactionsToBeRemoved))

ga = geneticA.GeneticAlgorithm(data, testCompStart, aminoAcids)        # initialise the GA with data
ga.population_size = 200                    # population size 
ga.generations = 3000		    # amount of generations


# define a fitness function
names = []
for rxn in (modelInit.reactions-reactionsToBeRemoved):
    names.append(rxn.id)

def fitness(individual, data, queue):
     modelGA = modelInit.copy()
     for i in range(len(individual.genes)):
         if individual.genes[i] == 0:
             if i < len(testCompStart):
                 testComp.append(testCompStart[i])
             else:
                 modelGA.reactions.get_by_id(names[i-len(testCompStart)]).lower_bound = 0
                 modelGA.reactions.get_by_id(names[i-len(testCompStart)]).upper_bound = 0
     score = evolveX(testComp,inh,modelGA,targets,dir,1)
     results = queue.get()
     results[individual.id] = score
     queue.put(results)
     

ga.fitness_function = fitness               # set the GA's fitness function
ga.run(testCompStart,names)                                    # run the GA
result = ga.best_individual()               # print the GA's best solution
bestGenes = ga.best_individualGenes()
print(result)
with open("solution_4.txt", 'w') as out:  #here solution name
	out.write(','.join(str(s) for s in result))
scores = ga.scoreslist
numberOfIterations = ga.generations
iterations = range(0,numberOfIterations)
plt.plot(iterations, scores, 'g', label='Scores')
plt.title('Scores during the GA')
plt.xlabel('Iteration number')
plt.ylabel('Score')
plt.savefig('scoreplot_4.png')  #here solution name

           
        
for i in range (len(bestGenes)):
    if bestGenes[i] == 0:
        if i < len(testCompStart):    
            print([i], testCompStart[i], bestGenes[i])
        else:
            print([i], [i-len(testCompStart)], names[i-len(testCompStart)], bestGenes[i])











 






