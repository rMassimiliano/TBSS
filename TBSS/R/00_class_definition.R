## This file pre-define the classes used in the package to avoid issues due to the fact that files are loaded in alphabetical order
#General TS class for tree scan
setClass(Class = "TS", slots = c(data = 'data.frame',
                                 tree ='data.frame',
                                 leaves = 'character',
                                 nodes = 'character',
                                 mapNodesLeaves = 'list',
				 nodesToTest ='character',
                                 nodeSS = 'data.frame',
                                 LRT ='numeric',
                                 B = 'numeric',
                                 LRT_H0 = 'matrix'))

setClass(Class = "unconditionalBernoulliTS", slots = c(p = 'numeric'), contains = 'TS')
setClass(Class = "unconditionalBernoulliTSVariableRR", contains = 'unconditionalBernoulliTS')


setClass(Class = "unconditionalPoissonTS",  contains = 'TS')

setClass(Class = "conditionalPoissonTS",  contains = 'unconditionalPoissonTS')
