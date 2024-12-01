#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import math
import sys
import z3
import time
import multiprocessing
import datetime
import glob


step = int(sys.argv[1])
mode_width = int(sys.argv[2])

maxAct = 2
maxRep = 2
threshold = 0.99
thresholdStep = 0.02
max_workers = int(sys.argv[3])

# binarized expression data
exppath = sys.argv[4]
exp_mat_binarized = pd.read_csv(exppath, delimiter='\t', index_col=0)

gene_names = exp_mat_binarized.columns.values
cell_names = exp_mat_binarized.index.values
print(exp_mat_binarized.shape)
print(exp_mat_binarized.iloc[0:5, 0:4])

# network dege data
edgepath = sys.argv[5]
network = pd.read_table(edgepath, sep="\t", header=0)

# read in list of genes to find logics
geneNames = gene_names

# Read in expression matrix
expression = exp_mat_binarized

# Read in path cells names as list
pathNames = cell_names

# into input/output names with desired mode width
mode_range = range(-int(mode_width / 2), int(mode_width / 2) + 1)
inputOutput = [[pathNames[i - step - j] for j in mode_range] + [pathNames[i]]
               for i in range(step + int(mode_width / 2), len(pathNames))]

inputOutputSelection = inputOutput

# How many genes are there
allGenes = list(expression.columns)
numGenes = len(allGenes)

# Define numbers corresponding to gates, genes and nothing variable
AND = 0
OR = 1
NOTHING = numGenes + 2

# from booleanRules_encodingFunctions import *
# Circuits are encoded as bitvectors
def makeCircuitVar(name):
    return z3.BitVec(name, 32)


def makeEnforcedVar(name):
    return z3.BitVec(name, 32)


# For each element in bitvector which gene/gate values is it allowed to take
def variableDomains(var, booleanAllowed, possibleInputs):
    if booleanAllowed == True:
        allowedValues = [0, 1] + possibleInputs
    else:
        allowedValues = possibleInputs

    return z3.Or([var == allowedValues[i] for i in range(len(allowedValues))])


# If a node = nothing then it is not allowed a parent as a gate
def parentsOfNothingArentGates(a, r):
    def f(c1, c2, p):
        return z3.Implies(z3.Or((c1 == NOTHING), (c2 == NOTHING)), z3.And(p != AND, p != OR))

    aParents = z3.And(
        z3.Implies(z3.Or(a[1] == NOTHING, a[2] == NOTHING), z3.And(a[0] != AND, a[0] != OR, a[0] != NOTHING)), \
        f(a[3], a[4], a[1]), \
        f(a[5], a[6], a[2]))

    rParents = z3.And(f(r[1], r[2], r[0]), \
                      f(r[3], r[4], r[1]), \
                      f(r[5], r[6], r[2]))

    return z3.And(aParents, rParents)


# If a node is a gene then it must have a gate as a parent
def parentsOfRestAreGates(a, r):
    def f(c1, c2, p):
        return z3.Implies(z3.Or((c1 != NOTHING), (c2 != NOTHING)), z3.Or(p == AND, p == OR))

    aParents = z3.And(f(a[1], a[2], a[0]), \
                      f(a[3], a[4], a[1]), \
                      f(a[5], a[6], a[2]))

    rParents = z3.And(f(r[1], r[2], r[0]), \
                      f(r[3], r[4], r[1]), \
                      f(r[5], r[6], r[2]))

    return z3.And(aParents, rParents)


# Can't have a gene more than once in the relation
def variablesDoNotAppearMoreThanOnce(symVars):
    def isVar(v):
        return z3.And(v != NOTHING, v != AND, v != OR)

    def notEqual(v, vars):
        return z3.And([v != i for i in vars if not z3.eq(v, i)])

    def doesNotAppearMoreThanOnce(v):
        return z3.Implies(isVar(v), notEqual(v, symVars))

    return z3.And([doesNotAppearMoreThanOnce(j) for j in symVars])


# Don't waste time doing things multiple times
def enforceSiblingLexigraphicalOrdering(v1, v2):
    return (v1 <= v2)


def enforceLexigraphicalOrderingBetweenBranches(p1, p2, c1, c2):
    return z3.Implies(p1 == p2, c1 <= c2)


def enforceLexigraphicalOrderingNaryGate(vars):
    return z3.Implies(vars[0] == vars[1], vars[2] <= vars[3])


# Store the activator and repressor variables in a list
activatorVars = ["a" + str(i) for i in range(7)]
repressorVars = ["r" + str(i) for i in range(7)]
circuitVars = activatorVars + repressorVars


# Depending on maximum number of inputs may want fewer nodes
def fixMaxInputs(v, max):
    if max == 0:
        return makeCircuitVar(v + "0") == NOTHING
    elif max == 1:
        return makeCircuitVar(v + "2") == NOTHING
    elif max == 2:
        return makeCircuitVar(v + "4") == NOTHING
    elif max == 3:
        return makeCircuitVar(v + "6") == NOTHING
    else:
        return True


def fixMaxActivators(max):
    return fixMaxInputs("a", max)


def fixMaxRepressors(max):
    return fixMaxInputs("r", max)


# This encodes the allowed update functions for a gene
def encodeUpdateFunction(gene, genes, maxActivators, maxRepressors, possAct, possRep):
    # Check all inputs are of right form
    assert (gene in genes and maxActivators > 0 and maxActivators <= 4 and maxRepressors >= 0 and maxRepressors <= 4), \
        "Incorrect arguments to encodeUpdateFunction"

    a = [makeCircuitVar("a%i" % i) for i in range(7)]
    r = [makeCircuitVar("r%i" % i) for i in range(7)]

    circuitEncoding = z3.And(variableDomains(a[0], True, possAct),
                             variableDomains(a[1], True, possAct + [NOTHING]),
                             variableDomains(a[2], True, possAct + [NOTHING]),
                             variableDomains(a[3], False, possAct + [NOTHING]),
                             variableDomains(a[4], False, possAct + [NOTHING]),
                             variableDomains(a[5], False, possAct + [NOTHING]),
                             variableDomains(a[6], False, possAct + [NOTHING]),
                             variableDomains(r[0], True, possRep + [NOTHING]),
                             variableDomains(r[1], True, possRep + [NOTHING]),
                             variableDomains(r[2], True, possRep + [NOTHING]),
                             variableDomains(r[3], False, possRep + [NOTHING]),
                             variableDomains(r[4], False, possRep + [NOTHING]),
                             variableDomains(r[5], False, possRep + [NOTHING]),
                             variableDomains(r[6], False, possRep + [NOTHING]), parentsOfNothingArentGates(a, r),
                             parentsOfRestAreGates(a, r), variablesDoNotAppearMoreThanOnce(a + r),
                             z3.And([enforceSiblingLexigraphicalOrdering(a[i], a[i + 1]) for i in [1, 3, 5]]),
                             z3.And([enforceSiblingLexigraphicalOrdering(r[i], r[i + 1]) for i in [1, 3, 5]]),
                             enforceLexigraphicalOrderingBetweenBranches(a[1], a[2], a[3], a[5]),
                             enforceLexigraphicalOrderingBetweenBranches(r[1], r[2], r[3], r[5]),
                             enforceLexigraphicalOrderingNaryGate(a), enforceLexigraphicalOrderingNaryGate(r),
                             fixMaxActivators(maxActivators), fixMaxRepressors(maxRepressors))

    return (circuitEncoding, a, r)


# Given inputs evaluates function
def evaluateUpdateFunction(aVars, rVars, geneValues, counter):
    i = counter

    intermediateValueVariablesA = [z3.Bool("va%i_%i" % (j, i)) for j in range(7)]
    intermediateValueVariablesR = [z3.Bool("vr%i_%i" % (j, i)) for j in range(7)]

    def andConstraints(symVars, variables, pi, c1i, c2i):
        return z3.Implies(symVars[pi] == z3.BitVecVal(AND, 32), variables[pi] == z3.And(variables[c1i], variables[c2i]))

    def orConstraints(symVars, variables, pi, c1i, c2i):
        return z3.Implies(symVars[pi] == z3.BitVecVal(OR, 32), variables[pi] == z3.Or(variables[c1i], variables[c2i]))

    def variableConstraints(symVars, intermediateVars):
        def f(symVar, i):
            return z3.And([z3.Implies(symVar == v, intermediateVars[i] == z3.BoolVal(geneValues[v - 2])) for v in
                           range(2, NOTHING)])

        return z3.And([f(symVars[i], i) for i in range(7)])

    circuitVal = z3.Bool("circuit_%i" % i)

    def circuitValue():
        noRepressors = rVars[0] == NOTHING

        return z3.If(noRepressors, intermediateValueVariablesA[0],
                     z3.And(intermediateValueVariablesA[0], z3.Not(intermediateValueVariablesR[0])))

    return (z3.And([variableConstraints(aVars, intermediateValueVariablesA),
                    variableConstraints(rVars, intermediateValueVariablesR),
                    andConstraints(aVars, intermediateValueVariablesA, 0, 1, 2),
                    andConstraints(aVars, intermediateValueVariablesA, 1, 3, 4),
                    andConstraints(aVars, intermediateValueVariablesA, 2, 5, 6),
                    andConstraints(rVars, intermediateValueVariablesR, 0, 1, 2),
                    andConstraints(rVars, intermediateValueVariablesR, 1, 3, 4),
                    andConstraints(rVars, intermediateValueVariablesR, 2, 5, 6),
                    orConstraints(aVars, intermediateValueVariablesA, 0, 1, 2),
                    orConstraints(aVars, intermediateValueVariablesA, 1, 3, 4),
                    orConstraints(aVars, intermediateValueVariablesA, 2, 5, 6),
                    orConstraints(rVars, intermediateValueVariablesR, 0, 1, 2),
                    orConstraints(rVars, intermediateValueVariablesR, 1, 3, 4),
                    orConstraints(rVars, intermediateValueVariablesR, 2, 5, 6),
                    circuitVal == circuitValue()]), circuitVal)


def circuitEvaluatesTo(gene, aVars, rVars, input, output, counter):
    outValue = output[gene - 2]
    inValues = input

    evaluationEncoding, circuitVal = evaluateUpdateFunction(aVars, rVars, inValues, counter)
    return (evaluationEncoding, circuitVal == z3.BoolVal(outValue))


# Define node class
class node:
    def __init__(self, pos_in, neg_in):
        self.p = pos_in
        self.n = neg_in


# Node list for network
nodeList = {}
nodeNames = list(set(network['to.gene']))

for n in nodeNames:
    allRelations = network[network['to.gene'] == n]
    posRelations = allRelations[allRelations['relation'] == 1]
    posGenes = list(posRelations['from.gene'])
    negRelations = allRelations[allRelations['relation'] == -1]
    negGenes = list(negRelations['from.gene'])
    nodeList[n] = node(posGenes, negGenes)

# Need to add some penalties in, remember she allowed self-activation
# Penalties
penAll = 0
penSelf = 0.005

# Dictionary to look up node names number equivalent
allNames = list(expression.columns)
nodeLookUp = {allNames[j]: j + 2 for j in range(len(allNames))}


# Add constraints to the solver
def constraintsBitVec(ctor, model, d):
    x = z3.BitVecVal(int(str(model[d])), 32)
    return ctor(str(d)) != x


def addConstraintsCircuitVar(solver, model, ds):
    constraints = z3.Or([constraintsBitVec(makeCircuitVar, model, d) for d in ds])
    return constraints


def modeCalc(lst):
    return max(set(lst), key=lst.count)


# Function which enforces minimum agreeing counts
def totalAgreement(inputOutputPair, gene, aVars, rVars, expressionValues, counter):
    outputName = inputOutputPair[-1]        # last element is output cell name
    inputName_list = inputOutputPair[:-1]   # input names are every element except the last one
    input_values = [expressionValues.loc[inputName,:] for inputName in inputName_list]
    input = [modeCalc([values[i] for values in input_values]) for i in range(len(input_values[0]))]
    output = expressionValues.loc[outputName, :]
    # counter += 1
    score = makeEnforcedVar("counter_%i" % counter)
    (encoding, match) = circuitEvaluatesTo(gene, aVars, rVars, input, output, counter)

    return (z3.And(encoding, z3.If(match, score == z3.BitVecVal(1, 32), score == z3.BitVecVal(0, 32))), score)


def agreesAtLeastNTimes(inputOutputList, expressionValues, gene, aVars, rVars, n):
    both = [totalAgreement(inputOutputList[p], gene, aVars, rVars, expressionValues, p) for p in
            range(len(inputOutputList))]
    encodings = [both[i][0] for i in range(len(both))]
    scoreValues = [both[i][1] for i in range(len(both))]

    return z3.And(z3.And(encodings), z3.Sum(scoreValues) >= n)


# Genes seems a bit pointless?
def findFunctions(solver, gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression,
                  agreementThreshold):
    expressionBool = expression == 1

    possAct = nodeList[gene].p
    possAct = [nodeLookUp[i] for i in possAct]
    possRep = nodeList[gene].n
    possRep = [nodeLookUp[i] for i in possRep]

    gene_str = gene
    gene = nodeLookUp[gene]
    genes = [nodeLookUp[genes[i]] for i in range(len(genes))]

    circuitEncoding, aVars, rVars = encodeUpdateFunction(gene, genes, maxActivators, maxRepressors, possAct, possRep)

    circuitVars = aVars + rVars

    # Choose number of agreement threshold
    agreementThreshold = int(agreementThreshold * len(inputOutput))

    solver.reset()
    allConstraints = z3.And(circuitEncoding, agreesAtLeastNTimes(inputOutput, expressionBool, gene, aVars, rVars, \
                                                                 agreementThreshold))
    solver.add(allConstraints)

    constraints = True

    possibleRules = []
    ruleNum = 0
    while str(solver.check()) == 'sat':
        m = solver.model()
        modelVariables = [m[i] for i in range(len(m))]
        circuitDecls = filter(lambda x: str(x) in [str(s) for s in circuitVars], modelVariables)
        enforceDecls = filter(lambda x: str(x)[0:7] == "counter", modelVariables)
        totalScore = sum([int(str(m[d])) for d in enforceDecls])

        newConstraints = addConstraintsCircuitVar(solver, m, circuitDecls)
        constraints = z3.And(constraints, newConstraints)

        solver.reset()
        solver.add(z3.And(allConstraints, constraints))

        possibleRules.append((m, totalScore))
        ruleNum += 1
        print('found rule %d' % ruleNum, 'for gene %s' % gene_str)

    if len(possibleRules) >= 500:
        newThreshold = max([s[1] for s in possibleRules])
        displayThreshold = float(newThreshold) / len(inputOutput)
        print('Finding too many rules. So now increase threshold to %f' % displayThreshold)

        solver.reset()
        allConstraints = z3.And(circuitEncoding, agreesAtLeastNTimes(inputOutput, expressionBool, gene, aVars, rVars, \
                                                                     newThreshold))
        solver.add(allConstraints)

        constraints = True
        possibleRules = []
        ruleNum = 0
        while str(solver.check()) == 'sat':
            m = solver.model()
            modelVariables = [m[i] for i in range(len(m))]
            circuitDecls = filter(lambda x: str(x) in [str(s) for s in circuitVars], modelVariables)
            enforceDecls = filter(lambda x: str(x)[0:7] == "counter", modelVariables)
            totalScore = sum([int(str(m[d])) for d in enforceDecls])

            newConstraints = addConstraintsCircuitVar(solver, m, circuitDecls)
            constraints = z3.And(constraints, newConstraints)

            solver.reset()
            solver.add(z3.And(allConstraints, constraints))

            possibleRules.append((m, totalScore))
            ruleNum += 1
            print('found rule %d' % ruleNum, 'for gene %s' % gene_str)
    return possibleRules, z3.And(allConstraints)


def findBestRuleForGene(gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression):
    s = z3.Solver()
    k = 0
    while 1==1:
        s.reset()
        agreement = threshold - thresholdStep * k
        print(datetime.datetime.now().strftime("%m%d - %H:%M:%S"),' Lowered threshold to %f' % agreement)
        rules, allConstraints = findFunctions(s, gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors,
                                              inputOutput, expression, agreement)
        k = k + 1
        if len(rules) > 0:
            print('Total %d rules' % len(rules), 'for gene %s' % gene); break
    convertedRules = [ruleConvert(r) for r in rules]
    geneInRules = []
    conRules = [x[0] for x in convertedRules]
    # for r in conRules:
    #     geneInRules.append(s[1] for s in r)
    #     import pdb
    #     # pdb.set_trace()
    #     if len(rules) > 0 and not all(g in r for r in geneInRules):
    #         break
    print('Found %d rules for gene satisfying threshold %f' % (len(rules), agreement))
    return rules, agreement, allConstraints


def checkAgreeLevel(solver, gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression,
                  agreementThreshold):
    expressionBool = expression == 1

    possAct = nodeList[gene].p
    possAct = [nodeLookUp[i] for i in possAct]
    possRep = nodeList[gene].n
    possRep = [nodeLookUp[i] for i in possRep]

    gene_str = gene
    gene = nodeLookUp[gene]
    genes = [nodeLookUp[genes[i]] for i in range(len(genes))]

    circuitEncoding, aVars, rVars = encodeUpdateFunction(gene, genes, maxActivators, maxRepressors, possAct, possRep)

    # Choose number of agreement threshold
    agreementThreshold = int(agreementThreshold * len(inputOutput))

    solver.reset()
    allConstraints = z3.And(circuitEncoding, agreesAtLeastNTimes(inputOutput, expressionBool, gene, aVars, rVars, \
                                                                 agreementThreshold))
    solver.add(allConstraints)

    solver_not_sat = not(str(solver.check()) == 'sat')
    return solver_not_sat


def findBestAgreeLevelForGene(gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression):
    s = z3.Solver()
    k = 0
    solver_not_sat = True
    while solver_not_sat:
        s.reset()
        agreement = threshold - thresholdStep * k
        print(datetime.datetime.now().strftime("%m%d - %H:%M:%S"),'Lowered threshold to %3f' % agreement)
        solver_not_sat = checkAgreeLevel(s, gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors,
                                              inputOutput, expression, agreement)
        k = k + 1
    return agreement




def ruleConvert(model):
    score = model[1] / len(inputOutput)
    modelRules = model[0]

    def int2gene(i):
        i = int(str(i))
        if i == AND:
            return 'and'
        elif i == OR:
            return 'or'
        else:
            return allNames[i - 2]

    namedVariables = ["a%i" % i for i in range(7)] + ["r%i" % i for i in range(7)]
    modelVariables = [modelRules[i] for i in range(len(modelRules))]
    usefulVariables = filter(lambda x: str(x) in namedVariables, modelVariables)
    modelConstraints = [(str(v), int2gene(modelRules[v])) for v in usefulVariables if
                        str(modelRules[v]) != str(NOTHING)]

    return modelConstraints, score


def evaluateRule(rule, input_list):
    # Remember to input boolean expression
    rules = {'a%d' % i: NOTHING for i in range(7)}
    rules.update({'r%d' % i: NOTHING for i in range(7)})

    convertBack = {'and': 0, 'or': 1}
    convertBack.update({gene: nodeLookUp[gene] for gene in allNames})

    for r in rule:
        rules[r[0]] = convertBack[r[1]]

    def getValue(v):
        return modeCalc([values[v-2] for values in input_list])

    # Find intermediate variables
    v = {'va%d' % i: getValue(rules['a%d' % i]) for i in range(7) if rules['a%d' % i] in range(2, NOTHING)}
    v.update({'vr%d' % i: getValue(rules['r%d' % i]) for i in range(7) if rules['r%d' % i] in range(2, NOTHING)})

    # Evaluate activators
    if rules['a1'] == 0:
        inter_a1 = v['va3'] and v['va4']
    elif rules['a1'] == 1:
        inter_a1 = v['va3'] or v['va4']
    elif rules['a1'] in range(2, NOTHING):
        inter_a1 = v['va1']
    else:
        inter_a1 = True

    if rules['a2'] == 0:
        inter_a2 = v['va5'] and v['va6']
    elif rules['a2'] == 1:
        inter_a2 = v['va5'] or v['va6']
    elif rules['a2'] in range(2, NOTHING):
        inter_a2 = v['va2']
    else:
        inter_a2 = True

    if rules['a0'] == 0:
        inter_a0 = inter_a1 and inter_a2
    elif rules['a0'] == 1:
        inter_a0 = inter_a1 or inter_a2
    else:
        inter_a0 = v['va0']

    # Evaluate repressors
    if rules['r0'] != NOTHING:
        if rules['r1'] == 0:
            inter_r1 = v['vr3'] and v['vr4']
        elif rules['r1'] == 1:
            inter_r1 = v['vr3'] or v['vr4']
        elif rules['r1'] in range(2, NOTHING):
            inter_r1 = v['vr1']
        else:
            inter_r1 = True

        if rules['r2'] == 0:
            inter_r2 = v['vr5'] and v['vr6']
        elif rules['r2'] == 1:
            inter_r2 = v['vr5'] or v['vr6']
        elif rules['r2'] in range(2, NOTHING):
            inter_r2 = v['vr2']
        else:
            inter_r2 = True

        if rules['r0'] == 0:
            inter_r0 = inter_r1 and inter_r2
        elif rules['r0'] == 1:
            inter_r0 = inter_r1 or inter_r2
        else:
            inter_r0 = v['vr0']

        return inter_a0 and not inter_r0
    else:
        return inter_a0



def logic_inference_main(g, DIR_NAME):
    print(datetime.datetime.now().strftime("%m%d - %H:%M:%S"), 'Trying to find rules for %s' % g)
    expressBool = expression == 1
    geneRules, agreement, solver = findBestRuleForGene(g, allNames, nodeLookUp, nodeList, maxAct, maxRep,
                                                       inputOutputSelection, expression)
    rulesForGene = [ruleConvert(r) for r in geneRules]
    print('Converted to readable format')

    def scoreForRule(rule):
        raw = 0
        for io in range(len(inputOutput)):
            outputName = inputOutput[io][-1]  # last element is output cell name
            inputName_list = inputOutput[io][:-1]  # input names are every element except the last one
            input_values = [expressBool.loc[inputName, :] for inputName in inputName_list]
            output = expressBool.loc[outputName, g]
            predictedOutput = evaluateRule(rule, input_values)
            if predictedOutput == output:
                raw += 1

        score = penalty = len(rule) * penAll
        if any(g in x for x in rule):
            penalty = penalty + (penSelf - penAll)
        score = raw - penalty

        # Might want to return more things?
        return (score, rule)

    scoreRules = [(scoreForRule(r[0]), 'z3 score %f' % r[1]) for r in rulesForGene]
    print('Found agreement level for each rule')
    bestRule = max(scoreRules, key=lambda x: x[0])
    maxValue = bestRule[0][0]
    allBest = []
    for rule in scoreRules:
        if rule[0][0] == maxValue:
            allBest.append(rule)

    print('The best rules are %s' % str(allBest))

    scoreRules = sorted(scoreRules, key=lambda x: x[0][0], reverse=True)

    print('Writing the rules to a file')

    if not os.path.exists(DIR_NAME):    # check result folder exist
        os.makedirs(DIR_NAME)           # if not, make one.
    Out_fname = DIR_NAME + '/%s_boolean_rules_%d.txt' % (g, step)
    f = open(Out_fname, 'w')
    f.write('Agreement level = %f \n' % agreement)
    f.write('The best rules were:\n')
    for item in allBest:
        f.write("%s\n" % str(item))
    f.write('Other rules found were:\n')
    for item in scoreRules:
        f.write("%s\n" % str(item))
    f.close()

    print('Found rules for %s' % g)


def agree_level_calculate(g, DIR_NAME):
    print(datetime.datetime.now().strftime("%m%d - %H:%M:%S"), 'Trying to find agreement level for %s' % g)
    agreement = findBestAgreeLevelForGene(g, allNames, nodeLookUp, nodeList, maxAct, maxRep,
                                                       inputOutputSelection, expression)
    print('Agreement level for gene %s is %f' %(g, agreement))
    if not os.path.exists(DIR_NAME):    # check result folder exist
        os.makedirs(DIR_NAME)           # if not, make one.
    f = open(os.path.join(DIR_NAME, '%s_%d_%.2f.txt' % (g, step, agreement)), 'w')
    f.write('Agreement level = %f \n' % agreement)
    f.close()


def write_whole_agree_level(dir_name):
    if not os.path.exists(dir_name):
        print("check path"); return
    os.chdir(dir_name)
    filenames = glob.glob("*.txt")
    agree_level_fname = "agree_level_whole.txt"
    f = open(agree_level_fname, 'w')
    for fname in filenames:
        gene_name = fname.split('_')[0]
        agree_level = float(fname.split('_')[2].split('.txt')[0])
        f.write('%s\t%.2f\n' %(gene_name, agree_level))
    f.close()
    os.chdir("../")     # into original dir


if __name__ == '__main__':
# def scns_main(input1, input2):
#     step = int(input1)
#     mode_width = int(input2)
    print("main started")
    main_start_time = time.time()


    print("started step %d width %d" % (step, mode_width))
    start_time = time.time()
    # DIR_NAME = "result_" + str(step) + "_" + str(mode_width)
    DIR_NAME = "logic_inf_results"
    def logic_inference_main_withDIR(g):
        logic_inference_main(g, DIR_NAME)
    def agree_level_calculate_withDIR(g):
        agree_level_calculate(g, DIR_NAME)

    pool = multiprocessing.Pool(processes=max_workers)
    pool.map(logic_inference_main_withDIR, gene_names)    # logic inference
    # pool.map(agree_level_calculate_withDIR, gene_names)     # only agreement level
    pool.close()
    pool.join()
    end_time = time.time()
    print("==========finished %d %d =========" % (step, mode_width))
    print("elapsed time : %f" % (end_time - start_time))
    # write_whole_agree_level(DIR_NAME)                       # not with logic inference

    print("==========finished everything =========")
    print("elapsed time : %f" % (time.time() - main_start_time))

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
